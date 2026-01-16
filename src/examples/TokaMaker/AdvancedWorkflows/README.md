# AdvancedWorkflows  
_Examples & reference implementations for advanced workflow scenarios utilizing TokaMaker_

## Overview  
This directory contains advanced example workflows built on top of the OpenFUSIONToolkit (OFT) framework within the `TokaMaker` example project. These workflows illustrate how to orchestrate multiple components, handle complex data flows, customize extension points, and integrate with external systems.

Use them as reference templates for your own custom pipelines, or as a base for extending the toolkit into production-ready scenarios.

## Prerequisites  
Before running any workflow in this folder, ensure the following:  
- You have cloned this repo and checked out the `main` (or appropriate) branch.  
- The OFT library version is compatible (specify version if applicable).  
- External systems or services referenced by workflows are accessible/configured (e.g., databases, message queues, cloud storage).  
- Environment variables or config files are set up (see next section).

## Folder Structure  
```text
AdvancedWorkflows/
│
├── Pulse Design/
│   ├── etc
│   ├── etc
│   └── README.md

├── TokaMax/
│   ├── etc
│   ├── etc
│   └── README.md
│
├── Poloidal Field Coil Optimization/
│   ├── etc
│   ├── etc
│   └── README.md
│
└── Isoflux Controller/
    ├── etc
    ├── etc
    └── README.md
