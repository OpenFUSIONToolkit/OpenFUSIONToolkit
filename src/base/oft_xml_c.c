/*-----------------------------------------------------------------------------
 * Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
 *
 * SPDX-License-Identifier: LGPL-3.0-only
 *------------------------------------------------------------------------------
 * Libxml2 DOM interface functions for the Open FUSION Toolkit
 *
 * Provides C functions callable from Fortran via iso_c_binding for
 * DOM-style access to XML files using the Libxml2 library.
 *----------------------------------------------------------------------------*/
#include <string.h>
#include <stdlib.h>
#include <libxml/parser.h>

/**
 * Parse an XML file from a given file path
 *
 * @param filepath Null-terminated path to XML file
 * @param doc_ptr  Output pointer to xmlDoc
 * @return 0 on success, nonzero on error
 */
int oft_xml_parse_file(const char* filepath, void** doc_ptr) {
    xmlDoc* doc;
    *doc_ptr = NULL;
    if (filepath == NULL) return 1;
    LIBXML_TEST_VERSION
    doc = xmlReadFile(filepath, NULL, 0);
    if (doc == NULL) return 2;
    *doc_ptr = (void*)doc;
    return 0;
}

/**
 * Retrieve a pointer to the root element of an XML document
 *
 * @param doc_ptr  Pointer to xmlDoc
 * @param root_ptr Output pointer to root xmlNode
 * @return 0 on success, nonzero on error
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

/**
 * Retrieve a pointer to the i-th xml node with a given name contained
 * within a specified parent node (1-based index)
 *
 * @param parent_ptr  Pointer to parent xmlNode
 * @param name        Null-terminated element name to search for
 * @param index       1-based index among matching children
 * @param element_ptr Output pointer to matching xmlNode
 * @return 0 on success, nonzero on error
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
                // printf("Found element: %s->%s\n", parent->name, cur->name);
                return 0;
            }
        }
    }
    return 2; /* not found */
}

/**
 * Retrieve pointers to all xml nodes with a given name contained within
 * a specified parent node
 *
 * @param parent_ptr   Pointer to parent xmlNode (cast to void*)
 * @param name         Null-terminated element name to search for
 * @param n            Output count of matching children
 * @param elements_ptr Output pointer to array of pointers to matching xmlNodes
 * @return 0 on success, nonzero on error
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

/**
 * Free a C array/object
 *
 * @param gen_ptr Pointer to be freed
 */
void oft_xml_free_ptr(void* gen_ptr) {
    if (gen_ptr != NULL) free(gen_ptr);
}

/**
 * Extract the string content from a given xml node.
 *
 * @param node_ptr    Pointer to xmlNode
 * @param content     Output character buffer (null-terminated)
 * @param content_len Size of content buffer
 * @return 0 on success, nonzero on error.
 */
int oft_xml_get_content(const void* node_ptr, char** content, int* content_len) {
    const xmlNode* node = (const xmlNode*)node_ptr;
    *content_len = 0;
    if (node == NULL) return 1;
    xmlChar* text = xmlNodeGetContent(node);
    if (text == NULL) {
        *content = NULL;
        return 0;
    }
    *content_len = (int)strlen((const char*)text);
    *content = (char*)malloc((*content_len + 1) * sizeof(char));
    strncpy(*content, (const char*)text, (size_t)*content_len);
    (*content)[*content_len] = '\0';
    xmlFree(text);
    *content_len += 1; /* include null terminator in length */
    return 0;
}

/**
 * Test if a given xml node has a specified attribute
 *
 * @param node_ptr    Pointer to xmlNode
 * @param attr_name   Null-terminated attribute name
 * @return 0 if attribute exists, nonzero otherwise
 */
int oft_xml_has_attribute(const void* node_ptr, const char* attr_name) {
    const xmlNode* node = (const xmlNode*)node_ptr;
    int len;
    if (node == NULL || attr_name == NULL) {
        return 1;
    }
    xmlChar* attr = xmlGetProp(node, (const xmlChar*)attr_name);
    if (attr == NULL) return 2;
    xmlFree(attr);
    return 0;
}

/**
 * Extract the string value of a given attribute on a given xml node
 *
 * @param node_ptr    Pointer to xmlNode
 * @param attr_name   Null-terminated attribute name
 * @param content     Output character buffer (null-terminated)
 * @param content_len Size of content buffer
 * @return 0 on success, nonzero on error.
 */
int oft_xml_get_attribute(const void* node_ptr, const char* attr_name,
                          char** content, int* content_len) {
    const xmlNode* node = (const xmlNode*)node_ptr;
    *content_len = 0;
    if (node == NULL || attr_name == NULL) {
        *content = NULL;
        return 1;
    }
    xmlChar* attr = xmlGetProp(node, (const xmlChar*)attr_name);
    if (attr == NULL) return 2;
    *content_len = (int)strlen((const char*)attr);
    *content = (char*)malloc((*content_len + 1) * sizeof(char));
    strncpy(*content, (const char*)attr, (size_t)*content_len);
    (*content)[*content_len] = '\0';
    xmlFree(attr);
    *content_len += 1; /* include null terminator in length */
    return 0;
}

/**
 * Free an XML document previously parsed with \ref oft_xml_load_file
 *
 * @param doc_ptr Pointer to xmlDoc
 */
void oft_xml_free_doc(void* doc_ptr) {
    if (doc_ptr != NULL)
        xmlFreeDoc((xmlDoc*)doc_ptr);
}