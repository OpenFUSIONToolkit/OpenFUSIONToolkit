/*-----------------------------------------------------------------------------
 * Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
 *
 * SPDX-License-Identifier: LGPL-3.0-only
 *------------------------------------------------------------------------------
 * LIBXML2 DOM interface functions for the Open FUSION Toolkit
 *
 * Provides C functions callable from Fortran via iso_c_binding for
 * DOM-style access to XML files using the LIBXML2 library.
 *----------------------------------------------------------------------------*/
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/*
 * Internal helper: trim leading/trailing whitespace from a token given by
 * pointer and length, then parse as a 32-bit integer.
 * Returns 0 on success, nonzero on error.
 */
static int parse_trim_int(const char* s, int len, int32_t* val) {
    char buf[128];
    int i, j, tlen;
    char* endptr;
    long v;
    if (len <= 0 || len >= 128) return 1;
    i = 0; while (i < len && isspace((unsigned char)s[i])) i++;
    j = len - 1; while (j > i && isspace((unsigned char)s[j])) j--;
    if (i > j) return 1;
    tlen = j - i + 1;
    memcpy(buf, s + i, (size_t)tlen);
    buf[tlen] = '\0';
    v = strtol(buf, &endptr, 10);
    if (*endptr != '\0') return 1;
    *val = (int32_t)v;
    return 0;
}

/*
 * Internal helper: trim leading/trailing whitespace from a token given by
 * pointer and length, then parse as a double.
 * Returns 0 on success, nonzero on error.
 */
static int parse_trim_real(const char* s, int len, double* val) {
    char buf[128];
    int i, j, tlen;
    char* endptr;
    if (len <= 0 || len >= 128) return 1;
    i = 0; while (i < len && isspace((unsigned char)s[i])) i++;
    j = len - 1; while (j > i && isspace((unsigned char)s[j])) j--;
    if (i > j) return 1;
    tlen = j - i + 1;
    memcpy(buf, s + i, (size_t)tlen);
    buf[tlen] = '\0';
    *val = strtod(buf, &endptr);
    if (*endptr != '\0') return 1;
    return 0;
}

/*
 * Internal helper: trim leading/trailing whitespace from a token given by
 * pointer and length, then parse as a logical (boolean) value.
 * Accepted strings (case-insensitive): T, TRUE, .TRUE., 1 (true);
 *                                       F, FALSE, .FALSE., 0 (false).
 * Returns 0 on success, nonzero on error.
 */
static int parse_trim_logical(const char* s, int len, int32_t* val) {
    char buf[16];
    int i, j, k, tlen;
    if (len <= 0 || len >= 16) return 1;
    i = 0; while (i < len && isspace((unsigned char)s[i])) i++;
    j = len - 1; while (j > i && isspace((unsigned char)s[j])) j--;
    if (i > j) return 1;
    tlen = j - i + 1;
    for (k = 0; k < tlen; k++) buf[k] = (char)toupper((unsigned char)s[i + k]);
    buf[tlen] = '\0';
    if (strcmp(buf, "T") == 0 || strcmp(buf, "TRUE") == 0 ||
            strcmp(buf, ".TRUE.") == 0 || strcmp(buf, "1") == 0) {
        *val = 1; return 0;
    }
    if (strcmp(buf, "F") == 0 || strcmp(buf, "FALSE") == 0 ||
            strcmp(buf, ".FALSE.") == 0 || strcmp(buf, "0") == 0) {
        *val = 0; return 0;
    }
    return 1;
}

/*
 * Parse a comma-and-newline-delimited string into a flat array of 32-bit
 * integers.  Commas separate columns within a row; newlines separate rows.
 * Empty/whitespace-only lines are skipped.
 *
 * str:    null-terminated input string
 * arr:    output flat array (caller-allocated, at least max_n elements)
 * max_n:  maximum number of elements to write
 * shape:  output 2-element array: shape[0]=nrows, shape[1]=ncols
 * Returns 0 on success, nonzero on error (-1 overflow, 2 parse error,
 *         3 inconsistent column count).
 */
int oft_xml_parse_int_array_c(const char* str, int32_t* arr, int max_n,
                               int32_t* shape) {
    const char* p;
    int nrows = 0, ncols_expected = -1, total = 0;
    shape[0] = 0; shape[1] = 0;
    if (str == NULL || arr == NULL || max_n <= 0) return 1;
    p = str;
    while (*p) {
        const char* row_start = p;
        int all_ws, ncols;
        const char* q;
        /* advance to end of row */
        while (*p && *p != '\n') p++;
        {
            const char* row_end = p;
            if (*p == '\n') p++;
            /* skip whitespace-only rows */
            all_ws = 1;
            {
                const char* c;
                for (c = row_start; c < row_end; c++) {
                    if (!isspace((unsigned char)*c)) { all_ws = 0; break; }
                }
            }
            if (all_ws) continue;
            /* parse comma-separated tokens in this row */
            ncols = 0;
            q = row_start;
            while (q <= row_end) {
                const char* tok_start = q;
                const char* tok_end;
                while (q < row_end && *q != ',') q++;
                tok_end = q;
                if (q < row_end) q++; /* skip comma */
                if (total >= max_n) return -1;
                if (parse_trim_int(tok_start, (int)(tok_end - tok_start),
                                   &arr[total]) != 0) return 2;
                total++;
                ncols++;
                if (tok_end >= row_end) break;
            }
            if (ncols_expected < 0) ncols_expected = ncols;
            else if (ncols != ncols_expected) return 3;
            nrows++;
        }
    }
    if (nrows == 0) return 1;
    shape[0] = (int32_t)nrows;
    shape[1] = (int32_t)ncols_expected;
    return 0;
}

/*
 * Parse a comma-and-newline-delimited string into a flat array of doubles.
 * See oft_xml_parse_int_array_c for format details.
 */
int oft_xml_parse_real_array_c(const char* str, double* arr, int max_n,
                                int32_t* shape) {
    const char* p;
    int nrows = 0, ncols_expected = -1, total = 0;
    shape[0] = 0; shape[1] = 0;
    if (str == NULL || arr == NULL || max_n <= 0) return 1;
    p = str;
    while (*p) {
        const char* row_start = p;
        int all_ws, ncols;
        const char* q;
        while (*p && *p != '\n') p++;
        {
            const char* row_end = p;
            if (*p == '\n') p++;
            all_ws = 1;
            {
                const char* c;
                for (c = row_start; c < row_end; c++) {
                    if (!isspace((unsigned char)*c)) { all_ws = 0; break; }
                }
            }
            if (all_ws) continue;
            ncols = 0;
            q = row_start;
            while (q <= row_end) {
                const char* tok_start = q;
                const char* tok_end;
                while (q < row_end && *q != ',') q++;
                tok_end = q;
                if (q < row_end) q++;
                if (total >= max_n) return -1;
                if (parse_trim_real(tok_start, (int)(tok_end - tok_start),
                                    &arr[total]) != 0) return 2;
                total++;
                ncols++;
                if (tok_end >= row_end) break;
            }
            if (ncols_expected < 0) ncols_expected = ncols;
            else if (ncols != ncols_expected) return 3;
            nrows++;
        }
    }
    if (nrows == 0) return 1;
    shape[0] = (int32_t)nrows;
    shape[1] = (int32_t)ncols_expected;
    return 0;
}

/*
 * Parse a comma-and-newline-delimited string into a flat array of logical
 * values encoded as int32_t (1=true, 0=false).
 * See oft_xml_parse_int_array_c for format details.
 */
int oft_xml_parse_logical_array_c(const char* str, int32_t* arr, int max_n,
                                   int32_t* shape) {
    const char* p;
    int nrows = 0, ncols_expected = -1, total = 0;
    shape[0] = 0; shape[1] = 0;
    if (str == NULL || arr == NULL || max_n <= 0) return 1;
    p = str;
    while (*p) {
        const char* row_start = p;
        int all_ws, ncols;
        const char* q;
        while (*p && *p != '\n') p++;
        {
            const char* row_end = p;
            if (*p == '\n') p++;
            all_ws = 1;
            {
                const char* c;
                for (c = row_start; c < row_end; c++) {
                    if (!isspace((unsigned char)*c)) { all_ws = 0; break; }
                }
            }
            if (all_ws) continue;
            ncols = 0;
            q = row_start;
            while (q <= row_end) {
                const char* tok_start = q;
                const char* tok_end;
                while (q < row_end && *q != ',') q++;
                tok_end = q;
                if (q < row_end) q++;
                if (total >= max_n) return -1;
                if (parse_trim_logical(tok_start, (int)(tok_end - tok_start),
                                       &arr[total]) != 0) return 2;
                total++;
                ncols++;
                if (tok_end >= row_end) break;
            }
            if (ncols_expected < 0) ncols_expected = ncols;
            else if (ncols != ncols_expected) return 3;
            nrows++;
        }
    }
    if (nrows == 0) return 1;
    shape[0] = (int32_t)nrows;
    shape[1] = (int32_t)ncols_expected;
    return 0;
}

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>

/*
 * Parse an XML file from a given file path.
 *
 * filepath: null-terminated path to XML file
 * doc_ptr:  output pointer to xmlDoc (cast to void*)
 * Returns 0 on success, nonzero on error.
 */
int oft_xml_load_file(const char* filepath, void** doc_ptr) {
    xmlDoc* doc;
    *doc_ptr = NULL;
    if (filepath == NULL) return 1;
    LIBXML_TEST_VERSION
    doc = xmlReadFile(filepath, NULL, 0);
    if (doc == NULL) return 2;
    *doc_ptr = (void*)doc;
    return 0;
}

/*
 * Retrieve a pointer to the root element of an XML document.
 *
 * doc_ptr:  pointer to xmlDoc (from oft_xml_load_file), cast to void*
 * root_ptr: output pointer to root xmlNode (cast to void*)
 * Returns 0 on success, nonzero on error.
 */
int oft_xml_get_root(const void* doc_ptr, void** root_ptr) {
    xmlNode* root;
    *root_ptr = NULL;
    if (doc_ptr == NULL) return 1;
    root = xmlDocGetRootElement((xmlDoc*)doc_ptr);
    if (root == NULL) return 2;
    *root_ptr = (void*)root;
    return 0;
}

/*
 * Retrieve a pointer to the I-th xml node with a given name contained
 * within a specified parent node (1-based index).
 *
 * parent_ptr:  pointer to parent xmlNode (cast to void*)
 * name:        null-terminated element name to search for
 * index:       1-based index among matching children
 * element_ptr: output pointer to matching xmlNode (cast to void*)
 * Returns 0 on success, nonzero on error.
 */
int oft_xml_get_element(const void* parent_ptr, const char* name, int index,
                        void** element_ptr) {
    const xmlNode* parent = (const xmlNode*)parent_ptr;
    xmlNode* cur;
    int count = 0;
    *element_ptr = NULL;
    if (parent == NULL || name == NULL || index <= 0) return 1;
    for (cur = parent->children; cur != NULL; cur = cur->next) {
        if (cur->type == XML_ELEMENT_NODE &&
                strcmp((const char*)cur->name, name) == 0) {
            count++;
            if (count == index) {
                *element_ptr = (void*)cur;
                return 0;
            }
        }
    }
    return 2; /* not found */
}

/*
 * Retrieve pointers to all xml nodes with a given name contained within
 * a specified parent node.
 *
 * parent_ptr:   pointer to parent xmlNode (cast to void*)
 * name:         null-terminated element name to search for
 * n:            output count of matching children
 * elements_ptr: output pointer to a heap-allocated array of void* pointers
 *               to matching xmlNodes; caller must free with oft_xml_free_elements
 * Returns 0 on success, nonzero on error.
 */
int oft_xml_get_elements(const void* parent_ptr, const char* name, int* n,
                         void** elements_ptr) {
    const xmlNode* parent = (const xmlNode*)parent_ptr;
    xmlNode* cur;
    void** arr;
    int count = 0;
    int i = 0;
    *n = 0;
    *elements_ptr = NULL;
    if (parent == NULL || name == NULL) return 1;
    /* Count matching children */
    for (cur = parent->children; cur != NULL; cur = cur->next) {
        if (cur->type == XML_ELEMENT_NODE &&
                strcmp((const char*)cur->name, name) == 0) {
            count++;
        }
    }
    if (count == 0) return 0;
    /* Allocate array of void* pointers */
    arr = (void**)malloc((size_t)count * sizeof(void*));
    if (arr == NULL) return 2;
    for (cur = parent->children; cur != NULL; cur = cur->next) {
        if (cur->type == XML_ELEMENT_NODE &&
                strcmp((const char*)cur->name, name) == 0) {
            arr[i++] = (void*)cur;
        }
    }
    *n = count;
    *elements_ptr = (void*)arr;
    return 0;
}

/*
 * Free the array of node pointers allocated by oft_xml_get_elements.
 *
 * elements_ptr: pointer to array previously returned by oft_xml_get_elements
 */
void oft_xml_free_elements(void* elements_ptr) {
    if (elements_ptr != NULL) free(elements_ptr);
}

/*
 * Extract the string content from a given xml node.
 *
 * node_ptr: pointer to xmlNode (cast to void*)
 * content:  output character buffer (not null-terminated by libxml2; we add '\0')
 * max_len:  maximum number of characters to write (including null terminator)
 * Returns 0 on success, nonzero on error.
 */
int oft_xml_get_content(const void* node_ptr, char* content, int max_len) {
    const xmlNode* node = (const xmlNode*)node_ptr;
    xmlChar* text;
    int len;
    if (node == NULL || content == NULL || max_len <= 0) return 1;
    text = xmlNodeGetContent(node);
    if (text == NULL) {
        content[0] = '\0';
        return 0;
    }
    len = (int)strlen((const char*)text);
    if (len >= max_len) len = max_len - 1;
    strncpy(content, (const char*)text, (size_t)len);
    content[len] = '\0';
    xmlFree(text);
    return 0;
}

/*
 * Extract the string value of a given attribute on a given xml node.
 *
 * node_ptr:  pointer to xmlNode (cast to void*)
 * attr_name: null-terminated attribute name
 * value:     output character buffer
 * max_len:   maximum number of characters to write (including null terminator)
 * Returns 0 on success, nonzero on error.
 */
int oft_xml_get_attribute(const void* node_ptr, const char* attr_name,
                          char* value, int max_len) {
    const xmlNode* node = (const xmlNode*)node_ptr;
    xmlChar* attr;
    int len;
    if (node == NULL || attr_name == NULL || value == NULL || max_len <= 0)
        return 1;
    attr = xmlGetProp(node, (const xmlChar*)attr_name);
    if (attr == NULL) return 2;
    len = (int)strlen((const char*)attr);
    if (len >= max_len) len = max_len - 1;
    strncpy(value, (const char*)attr, (size_t)len);
    value[len] = '\0';
    xmlFree(attr);
    return 0;
}

/*
 * Free an XML document previously parsed with oft_xml_load_file.
 *
 * doc_ptr: pointer to xmlDoc (from oft_xml_load_file), cast to void*
 */
void oft_xml_free_doc(void* doc_ptr) {
    if (doc_ptr != NULL)
        xmlFreeDoc((xmlDoc*)doc_ptr);
}

#endif /* HAVE_LIBXML2 */
