#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-3.0-only
"""Install the OFT agent skills and subagents for Claude Code, Codex, or Gemini CLI.

Skills (.agents/skills/<name>/, Agent Skills format) are copied as-is (or symlinked
with --link). Subagents are kept in Claude Code format (.claude/agents/<name>.md) and
translated to each tool's format, with models mapped through model_map.json.

Usage:
  python .agents/install.py --target anthropic --scope user
  python .agents/install.py --target openai google --scope project --project-dir ~/my_study
  python .agents/install.py --list
  python .agents/install.py --target openai --scope user --uninstall

Destinations:
  anthropic (Claude Code)  skills: .claude/skills/   agents: .claude/agents/<name>.md
  openai    (Codex)        skills: .agents/skills/   agents: .codex/agents/<name>.toml
  google    (Gemini CLI)   skills: .agents/skills/   agents: .gemini/agents/<name>.md
under the home directory (--scope user) or the project directory (--scope project).
"""
import argparse
import json
import os
import re
import shutil
import sys

SRC_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKILLS_SRC = os.path.join(SRC_ROOT, '.agents', 'skills')
AGENTS_SRC = os.path.join(SRC_ROOT, '.claude', 'agents')
MODEL_MAP_FILE = os.path.join(SRC_ROOT, '.agents', 'model_map.json')

# Skills installed when --skills is not given: the ones for running the codes. The
# others (oft-build, oft-test) are for developing OFT and live in the repository.
DEFAULT_SKILLS = ['tokamaker-equilibrium', 'thincurr-model', 'thincurr-run']

TARGET_ALIASES = {'anthropic': 'anthropic', 'claude': 'anthropic',
                  'openai': 'openai', 'codex': 'openai',
                  'google': 'google', 'gemini': 'google'}
TARGETS = {
    'anthropic': {'label': 'Claude Code', 'skills': ('.claude', 'skills'), 'agents': ('.claude', 'agents'), 'ext': '.md'},
    'openai': {'label': 'Codex', 'skills': ('.agents', 'skills'), 'agents': ('.codex', 'agents'), 'ext': '.toml'},
    'google': {'label': 'Gemini CLI', 'skills': ('.agents', 'skills'), 'agents': ('.gemini', 'agents'), 'ext': '.md'},
}

# Claude Code tool names -> Gemini CLI tool names
GEMINI_TOOLS = {'Read': 'read_file', 'Write': 'write_file', 'Edit': 'replace', 'MultiEdit': 'replace',
                'Bash': 'run_shell_command', 'Glob': 'glob', 'Grep': 'grep_search', 'LS': 'list_directory',
                'WebFetch': 'web_fetch', 'WebSearch': 'google_web_search', 'TodoWrite': 'write_todos'}
CLAUDE_FAMILIES = ('fable', 'opus', 'sonnet', 'haiku')
COPY_IGNORE = shutil.ignore_patterns('__pycache__', '*.pyc', '.ruff_cache', '.DS_Store')


# ---------------------------------------------------------------------------
# Reading the sources
# ---------------------------------------------------------------------------
def parse_frontmatter(text, filename):
    """Split a Markdown file with YAML frontmatter into (fields, body).

    Handles the subset used by the agent files: "key: scalar" lines and "key:"
    followed by "  - item" lines.
    """
    lines = text.splitlines()
    if not lines or lines[0].strip() != '---':
        raise ValueError('{0}: no YAML frontmatter'.format(filename))
    try:
        end = lines.index('---', 1)
    except ValueError:
        raise ValueError('{0}: unterminated YAML frontmatter'.format(filename))
    fields = {}
    key = None
    for line in lines[1:end]:
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        if line.lstrip().startswith('- ') and key is not None:
            if not isinstance(fields[key], list):
                fields[key] = []
            fields[key].append(_scalar(line.lstrip()[2:]))
            continue
        if ':' not in line or line[0] in ' \t':
            raise ValueError('{0}: unsupported frontmatter line "{1}"'.format(filename, line))
        key, value = line.split(':', 1)
        key = key.strip()
        fields[key] = _scalar(value) if value.strip() else ''
    body = '\n'.join(lines[end+1:]).strip('\n') + '\n'
    return fields, body


def _scalar(value):
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in '"\'':
        return value[1:-1]
    return value


def read_agents():
    """Read the source agent definitions, keyed by name."""
    agents = {}
    if not os.path.isdir(AGENTS_SRC):
        return agents
    for filename in sorted(os.listdir(AGENTS_SRC)):
        if not filename.endswith('.md'):
            continue
        path = os.path.join(AGENTS_SRC, filename)
        with open(path, 'r') as fid:
            fields, body = parse_frontmatter(fid.read(), path)
        for required in ('name', 'description'):
            if not fields.get(required):
                raise ValueError('{0}: missing "{1}"'.format(path, required))
        for list_key in ('tools', 'skills'):
            if isinstance(fields.get(list_key), str):
                fields[list_key] = [item.strip() for item in fields[list_key].split(',') if item.strip()]
        agents[fields['name']] = {'fields': fields, 'body': body, 'path': path}
    return agents


def available_skills():
    """Names of the skills in the repository (directories holding a SKILL.md)."""
    return sorted(name for name in os.listdir(SKILLS_SRC)
                  if os.path.isfile(os.path.join(SKILLS_SRC, name, 'SKILL.md')))


def load_model_map(extra_file):
    """Read model_map.json, then overlay the entries from --model-map, if given."""
    with open(MODEL_MAP_FILE, 'r') as fid:
        model_map = json.load(fid)
    if extra_file is not None:
        with open(extra_file, 'r') as fid:
            extra = json.load(fid)
        for target, entries in extra.items():
            if target.startswith('_'):
                continue
            target = TARGET_ALIASES.get(target, target)
            model_map.setdefault(target, {}).update(entries)
    return model_map


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------
def map_model(model, target, model_map, warnings, agent_name):
    """Return the target's model for a source model value, or None to inherit."""
    if not model or model == 'inherit':
        return None
    entries = model_map.get(target, {})
    if model in entries:
        return entries[model]
    for family in CLAUDE_FAMILIES:
        if family in model.lower() and family in entries:
            return entries[family]
    if target == 'anthropic':
        return model
    warnings.append('{0}: no {1} mapping for model "{2}"; the agent will inherit the main model '
                    '(add it to model_map.json or --model-map)'.format(agent_name, target, model))
    return None


def agent_name_for(name, target):
    """Codex agent names follow its documentation's convention (underscores)."""
    return name.replace('-', '_') if target == 'openai' else name


def rename_agents(text, target, names):
    """Rewrite references to the agents in a text to the target's agent names."""
    for name in names:
        pattern = r'(?<![\w-]){0}(?![\w-])'.format(re.escape(name))
        text = re.sub(pattern, agent_name_for(name, target), text)
    return text


def skill_preamble(skills, target):
    """Instructions that load the agent's skills where they cannot be preloaded."""
    if not skills:
        return ''
    names = ', '.join('`{0}`'.format(skill) for skill in skills)
    if target == 'google':
        return ('Before anything else, call `activate_skill` for the {0} skill{1} and follow '
                'it.\n\n'.format(names, 's' if len(skills) > 1 else ''))
    return ('Before anything else, load the {0} skill{1} (`${2}`): read its `SKILL.md` and '
            'follow it.\n\n'.format(names, 's' if len(skills) > 1 else '', skills[0]))


def yaml_scalar(value):
    """A YAML scalar: plain for simple tokens, else double-quoted (JSON syntax is valid YAML)."""
    text = str(value)
    if text and all(c.isalnum() or c in '-_.' for c in text):
        return text
    return json.dumps(text)


def translate_anthropic(agent, model_map, warnings, names):
    fields = dict(agent['fields'])
    model = map_model(fields.get('model'), 'anthropic', model_map, warnings, fields['name'])
    if model is not None:
        fields['model'] = model
    lines = ['---']
    for key, value in fields.items():
        if isinstance(value, list):
            lines.append('{0}:'.format(key))
            lines.extend('  - {0}'.format(yaml_scalar(item)) for item in value)
        else:
            lines.append('{0}: {1}'.format(key, yaml_scalar(value)))
    lines.append('---')
    return '\n'.join(lines) + '\n\n' + agent['body']


def translate_google(agent, model_map, warnings, names):
    fields = agent['fields']
    name = fields['name']
    lines = ['---', 'name: {0}'.format(yaml_scalar(name)), 'description: {0}'.format(yaml_scalar(fields['description']))]
    if fields.get('tools'):
        tools = []
        for tool in fields['tools']:
            if tool in GEMINI_TOOLS:
                if GEMINI_TOOLS[tool] not in tools:
                    tools.append(GEMINI_TOOLS[tool])
            else:
                warnings.append('{0}: no Gemini CLI equivalent for tool "{1}"; dropped'.format(name, tool))
        if fields.get('skills'):
            tools.append('activate_skill')
        lines.append('tools:')
        lines.extend('  - {0}'.format(tool) for tool in tools)
    model = map_model(fields.get('model'), 'google', model_map, warnings, name)
    if model is not None:
        lines.append('model: {0}'.format(yaml_scalar(model)))
    if fields.get('maxTurns'):
        lines.append('max_turns: {0}'.format(int(fields['maxTurns'])))
    lines.append('---')
    body = skill_preamble(fields.get('skills'), 'google') + rename_agents(agent['body'], 'google', names)
    return '\n'.join(lines) + '\n\n' + body


def toml_string(value):
    return json.dumps(str(value), ensure_ascii=False)


def translate_openai(agent, model_map, warnings, names):
    fields = agent['fields']
    name = fields['name']
    body = skill_preamble(fields.get('skills'), 'openai') + rename_agents(agent['body'], 'openai', names)
    tools = fields.get('tools')
    if tools and not ({'Write', 'Edit', 'MultiEdit'} & set(tools)):
        body += '\nDo not create or modify files: report what should change instead.\n'
    lines = ['name = {0}'.format(toml_string(agent_name_for(name, 'openai'))),
             'description = {0}'.format(toml_string(rename_agents(fields['description'], 'openai', names)))]
    model = map_model(fields.get('model'), 'openai', model_map, warnings, name)
    if model is not None:
        lines.append('model = {0}'.format(toml_string(model)))
    if fields.get('effort'):
        lines.append('model_reasoning_effort = {0}'.format(toml_string(fields['effort'])))
    if "'''" in body:
        lines.append('developer_instructions = {0}'.format(toml_string(body)))
    else:
        lines.append("developer_instructions = '''\n{0}'''".format(body))
    if tools and 'Bash' not in tools:
        lines += ['', '[features]', 'shell_tool = false']
    return '\n'.join(lines) + '\n'


TRANSLATORS = {'anthropic': translate_anthropic, 'openai': translate_openai, 'google': translate_google}


# ---------------------------------------------------------------------------
# Installing
# ---------------------------------------------------------------------------
def same_path(path_a, path_b):
    return os.path.realpath(path_a) == os.path.realpath(path_b)


def remove_path(path):
    """Remove a file, a symlink (never its target), or a directory tree."""
    if os.path.islink(path) or os.path.isfile(path):
        os.unlink(path)
    else:
        shutil.rmtree(path)


def plan(args, model_map, agents):
    """List the (action, destination, source/content, description) steps to perform."""
    base = os.path.expanduser('~') if args.scope == 'user' else os.path.abspath(args.project_dir)
    steps = []
    warnings = []
    seen = set()
    for target in args.target:
        info = TARGETS[target]
        skill_dir = os.path.join(base, *info['skills'])
        for skill in args.skills:
            dest = os.path.join(skill_dir, skill)
            if dest in seen:
                continue
            seen.add(dest)
            steps.append(('skill', dest, os.path.join(SKILLS_SRC, skill), '{0} skill'.format(skill)))
        agent_dir = os.path.join(base, *info['agents'])
        for name in args.agents:
            agent = agents[name]
            missing = [skill for skill in agent['fields'].get('skills', []) if skill not in args.skills]
            if missing:
                warnings.append('{0}: skipped, it needs skill(s) not being installed: {1}'.format(name, ', '.join(missing)))
                continue
            dest = os.path.join(agent_dir, agent_name_for(name, target) + info['ext'])
            content = TRANSLATORS[target](agent, model_map, warnings, list(agents))
            steps.append(('agent', dest, (agent['path'], content), '{0} agent ({1})'.format(name, info['label'])))
    return steps, warnings


def apply_step(step, args):
    """Perform one install/uninstall step. Returns a one-line status."""
    kind, dest, source, label = step
    src_path = source if kind == 'skill' else source[0]
    exists = os.path.lexists(dest)
    is_link = os.path.islink(dest)
    # Never touch the repository's own copy (eg. a project install into the OFT repository).
    # A symlink is not a copy: unlinking it leaves its target alone.
    if exists and (not is_link) and same_path(dest, src_path):
        return 'in place   {0}: {1} (the repository copy)'.format(label, dest)
    if exists and is_link and args.link and (not args.uninstall) and same_path(dest, src_path):
        return 'in place   {0}: {1} (already linked)'.format(label, dest)
    if args.uninstall:
        if not exists:
            return 'absent     {0}: {1}'.format(label, dest)
        if not args.dry_run:
            remove_path(dest)
        return 'removed    {0}: {1}'.format(label, dest)
    if exists and not args.force:
        return 'exists     {0}: {1} (use --force to replace)'.format(label, dest)
    if not args.dry_run:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if exists:
            remove_path(dest)
        if kind == 'skill' and args.link:
            os.symlink(src_path, dest)
        elif kind == 'skill':
            shutil.copytree(src_path, dest, ignore=COPY_IGNORE)
        else:
            with open(dest, 'w') as fid:
                fid.write(source[1])
    verb = 'linked' if (kind == 'skill' and args.link) else ('replaced' if exists else 'installed')
    return '{0:10} {1}: {2}'.format(verb, label, dest)


def ask(prompt, choices):
    """Prompt for one of several choices (interactive use only)."""
    while True:
        answer = input('{0} [{1}]: '.format(prompt, '/'.join(choices))).strip().lower()
        if answer in choices:
            return answer


def main():
    parser = argparse.ArgumentParser(description='Install the OFT agent skills and subagents for Claude Code (anthropic), '
                                     'Codex (openai), or Gemini CLI (google).')
    parser.add_argument('--target', nargs='+', choices=sorted(TARGET_ALIASES), default=None,
                        help='Tool(s) to install for: anthropic/claude, openai/codex, google/gemini')
    parser.add_argument('--scope', choices=('user', 'project'), default=None,
                        help='Install under your home directory (user) or a project directory (project)')
    parser.add_argument('--project-dir', default='.', help='Project directory for --scope project (default: current directory)')
    parser.add_argument('--skills', nargs='+', default=None,
                        help='Skills to install (default: {0}); "all" for every skill'.format(' '.join(DEFAULT_SKILLS)))
    parser.add_argument('--agents', nargs='*', default=None,
                        help='Agents to install (default: every agent whose skills are installed); give none to skip agents')
    parser.add_argument('--model-map', default=None, help='JSON file with model-map entries overriding model_map.json')
    parser.add_argument('--link', action='store_true', help='Symlink skills to this repository instead of copying them')
    parser.add_argument('--force', action='store_true', help='Replace existing skills/agents with the same names')
    parser.add_argument('--uninstall', action='store_true', help='Remove the selected skills/agents from the destination')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without changing anything')
    parser.add_argument('--list', action='store_true', help='List the available skills and agents, then exit')
    args = parser.parse_args()

    agents = read_agents()
    skills = available_skills()
    if args.list:
        print('Skills ({0}):'.format(SKILLS_SRC))
        for skill in skills:
            print('  {0}{1}'.format(skill, '  (default)' if skill in DEFAULT_SKILLS else ''))
        print('Agents ({0}):'.format(AGENTS_SRC))
        for name, agent in agents.items():
            print('  {0}  (skills: {1}; model: {2})'.format(name, ', '.join(agent['fields'].get('skills', [])) or '-',
                                                        agent['fields'].get('model', 'inherit')))
        return 0

    interactive = sys.stdin.isatty()
    if args.target is None:
        if not interactive:
            parser.error('--target is required')
        args.target = [ask('Install for which tool', ['anthropic', 'openai', 'google'])]
    if args.scope is None:
        if not interactive:
            parser.error('--scope is required')
        args.scope = ask('Install for your user account (home directory) or a project', ['user', 'project'])
        if args.scope == 'project' and args.project_dir == '.':
            args.project_dir = input('Project directory [{0}]: '.format(os.getcwd())).strip() or '.'
    args.target = list(dict.fromkeys(TARGET_ALIASES[target] for target in args.target))
    if args.scope == 'project' and not os.path.isdir(args.project_dir):
        parser.error('project directory "{0}" does not exist'.format(args.project_dir))

    if args.skills is None:
        args.skills = list(DEFAULT_SKILLS)
    elif args.skills == ['all']:
        args.skills = skills
    unknown = [skill for skill in args.skills if skill not in skills]
    if unknown:
        parser.error('unknown skill(s): {0} (available: {1})'.format(', '.join(unknown), ', '.join(skills)))
    if args.agents is None:
        args.agents = [name for name, agent in agents.items()
                       if all(skill in args.skills for skill in agent['fields'].get('skills', []))]
    unknown = [name for name in args.agents if name not in agents]
    if unknown:
        parser.error('unknown agent(s): {0} (available: {1})'.format(', '.join(unknown), ', '.join(agents)))

    model_map = load_model_map(args.model_map)
    steps, warnings = plan(args, model_map, agents)
    if args.dry_run:
        print('Dry run: nothing is changed')
    for step in steps:
        print(apply_step(step, args))
    for warning in warnings:
        print('WARNING: ' + warning)

    if not args.uninstall:
        print('\nNext steps:')
        print('- The skills need a working OFT Python install: set OFT_ROOTPATH to its prefix, or put its')
        print('  python/ directory on PYTHONPATH (see each skill\'s "Before anything else").')
        if 'anthropic' in args.target:
            print('- Claude Code: restart sessions to pick up new skills and agents.')
        if 'openai' in args.target:
            print('- Codex: restart to pick up new skills and agents; ask for an agent by name (eg. "use thincurr_runner").')
        if 'google' in args.target:
            print('- Gemini CLI: restart; project skills and agents load only in trusted folders, and each project')
            print('  agent must be acknowledged when first seen.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
