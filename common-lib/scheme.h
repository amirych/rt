#ifndef SCHEME_H
#define SCHEME_H

#include<string>
#include<list>
#include<iostream>
#include<fstream>
#include<stdint.h>

#include<mxml.h>

#include "surface.h"
#include "beam.h"
#include "tabulated_function.h"

using namespace std;


#define ATTR_NAME_MAXLEN 1024


class Scheme
{
public:
    enum SchemeError {FileOpenFailure,ParseFailure,NoClassAttr,BadAttrValue,NanAttrValue,MemAllocFailure,
                      NoSchemeDesc,NoBeamDesc,NoSurfaceDesc,NoRange,NoRangeFile,BadRange,
                      EmptyScheme,BadStartIndex,FileSaveFailure};

    Scheme();
    Scheme(const string &filename);
    ~Scheme();

    // if QE_flag is true then read only scheme attributes (not all!) needed for QE computation!
    // the call the function deletes existing internal scheme description!
    void load_from_file(const string &filename, bool QE_flag = false);
    void destroy();
    void run(list<string> &resultfile_list, ostream &log_file = std::cout, const size_t start_surface = 0);

    size_t GetNumberOfSurfaces() const;

    void ComputeQE(size_t Nlambda, char* QE_filename);

    // Scheme class exceptions
    class bad_scheme: public exception
    {
    public:
        bad_scheme(SchemeError err);
        SchemeError scheme_error() const;
    private:
        SchemeError Error;
    };
private:
    Beam *beam;
    list<Surface*> scheme;
    list<vector<real_t> > wavelenth_ranges;

    string SchemeFileName;

    string SchemeName;
    string ResultFileRoot;
    string LogFileName;

    mxml_node_t *scheme_tree, *current_node;
    mxml_index_t *surface_index, *beam_index, *desc_index;

    const char* get_attr_value(string attr_name);
    string get_string_attr(string attr_name, bool toupper_case = true);
    vector<real_t> get_numeric_attr(string attr_name, bool strict = false);

    list<vector<real_t> > get_range(string range_attr_name, uint8_t number_of_grating);

    vector<real_t> total_QE_curve;
};

#endif // SCHEME_H
