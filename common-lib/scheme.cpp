#include "scheme.h"
#include "ascii_file.h"
#include "rt_engine_func.h"
#include "rt_engine_errors.h"
#include "version.h"
#include "num_defs.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mxml.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include <algorithm>

#include <stdint.h>

#ifdef USING_MSVC
    #include <sys/types.h>
    #include <sys/stat.h>
#endif

/* Auxiliary functions */

// function deletes leading and trailing whitesaces
static void trim_spaces(std::string& str, const std::string& whitespace = " \t")
{
    std::size_t strBegin = str.find_first_not_of(whitespace);

    if (strBegin == std::string::npos) {
        str.clear();
        return; // empty string
    }

    std::size_t strEnd = str.find_last_not_of(whitespace);
    std::size_t strRange = strEnd - strBegin + 1;

    str.assign(str, strBegin, strRange);
}

// function transform a string to upper case
static void toupcase(string &str) {
    for ( size_t i = 0; i < str.size(); ++i ) str[i] = toupper(str[i]);
}



         /********************************
         *                               *
         *  Scheme class implementation  *
         *                               *
         ********************************/


        /*  Scheme class exceptions implementation */

Scheme::bad_scheme::bad_scheme(SchemeError err): exception(), Error(err)
{

}

Scheme::SchemeError Scheme::bad_scheme::scheme_error() const
{
    return Error;
}


        /*  Scheme class constructors and destructor  */

Scheme::Scheme(): scheme(list<Surface*>()), beam(NULL), wavelenth_ranges(list<vector<real_t> >()),
    SchemeFileName(""),SchemeName("Unknown"), ResultFileRoot(""), LogFileName(""), total_QE_curve(vector<real_t>())
{
}


Scheme::Scheme(const string &filename): scheme(list<Surface*>()), beam(NULL), wavelenth_ranges(list<vector<real_t> >()),
    SchemeFileName(filename),SchemeName("Unknown"), ResultFileRoot(""), LogFileName(""), total_QE_curve(vector<real_t>())
{
    try {
        load_from_file(filename);
    } catch (bad_scheme &ex) {
        SchemeError err = ex.scheme_error();
        switch ( err ) {
            case Scheme::FileOpenFailure: {
                cerr << "Can not open input scheme file!!!\n";
                break;
            }
            case Scheme::ParseFailure: {
                cerr << "Can not parse input scheme desciption!!!\n";
                break;
            }
            case Scheme::NoClassAttr: {
                cerr << "There is no surface class description!!!\n";
                break;
            }
            case Scheme::BadAttrValue: {
                cerr << "Invalid surface attribute in the input scheme description!!!\n";
                break;
            }
            case Scheme::MemAllocFailure: {
                cerr << "Memory allocation failure!!!\n";
                break;
            }
            default: cerr << "Unknown error!!!\n";
        }
////        // if error acqured delete XML tree and index
//        mxmlIndexDelete(surface_index);
//        mxmlDelete(scheme_tree);
        // destroy scheme if exists
        destroy();
    } catch (...) {
        throw;
    }
}


Scheme::~Scheme()
{
    destroy();
}


        /*  Scheme class private methods  */

//const char*  Scheme::get_attr_value(mxml_node_t *node,string attr_name)
const char*  Scheme::get_attr_value(string attr_name)
{
    const char* attr_value;

//    mxml_node_t *attr = mxmlFindPath(node,attr_name.c_str());
    mxml_node_t *attr = mxmlFindPath(current_node,attr_name.c_str());

    if ( attr == NULL ) {
//        cerr << "ATTR: " << attr_name.c_str() << endl;
        throw bad_scheme(Scheme::ParseFailure);
    }

    attr_value = mxmlGetOpaque(attr);
    if ( attr_value == NULL ) {
        throw bad_scheme(Scheme::BadAttrValue);
    }

    return attr_value;
}


//string Scheme::get_string_attr(mxml_node_t *node, string attr_name)
string Scheme::get_string_attr(string attr_name, bool toupper_case)
{
    const char* attr_value;
    string str_value;

//    attr_value = get_attr_value(node, attr_name);
    try {
        attr_value = get_attr_value(attr_name);
    } catch (Scheme::bad_scheme &ex) { // attribute in the input file was not found
        throw;
    }

    str_value = attr_value;

    trim_spaces(str_value);
    if ( str_value.empty() ) {
        throw bad_scheme(Scheme::BadAttrValue);
    }

    if ( toupper_case ) toupcase(str_value);

    return str_value;
}


//vector<real_t> Scheme::get_numeric_attr(mxml_node_t *node, string attr_name)
vector<real_t> Scheme::get_numeric_attr(string attr_name, bool strict)
{
    const char* attr_value;
    char *p;
    string str_value;
    vector<real_t> nums;
    real_t val;

//    str_value = get_string_attr(node, attr_name);
    try {
        str_value = get_string_attr(attr_name);
    } catch (Scheme::bad_scheme &ex) {
        throw;
    }

    attr_value = str_value.c_str();

    // try to convert to doubles
    for (;;) {
        val = strtod(attr_value,&p);
        if ( p!= attr_value ) nums.push_back(val); else break;
        attr_value = p;
    }

    if ( strict ) { // check whether string did contain strictly a number or numeric vector
        if ( attr_value[0] != '\0') {
            throw bad_scheme(Scheme::NanAttrValue);
        }
    }

    if ( nums.empty() ) throw bad_scheme(Scheme::BadAttrValue);

    return nums;
}


list<vector<real_t> > Scheme::get_range(string range_attr_name, uint8_t number_of_grating)
{
    vector<real_t> val;
    string str_val;
    list<vector<real_t> > ranges;
    AsciiFile file;

    try {
        val = get_numeric_attr(range_attr_name,true); // is numeric value?
        if ( (val.size() < 2) || (val[0] >= val[1] ) || (val[1] <= 0.0) ) {
            throw Scheme::bad_scheme(Scheme::BadRange);
        }
        ranges.push_back(val);
    } catch (Scheme::bad_scheme &ex) {
        if ( ex.scheme_error() != Scheme::NanAttrValue ) throw; // attribute did not found or it is bad

        str_val = get_string_attr(range_attr_name,false); // a filename is given

        try { // read file
            file.open(str_val.c_str());

            if ( file.fail() ) {
                throw Scheme::bad_scheme(Scheme::NoRangeFile);
            }

            AsciiFile::AsciiFileFlag line_flag;

            while ( (line_flag = file.ReadLine(2+number_of_grating,val)) != AsciiFile::Eof ) {
                if ( line_flag == AsciiFile::InvalidData ) throw Scheme::bad_scheme(Scheme::BadRange);
                if ( line_flag == AsciiFile::DataString ) {
                    if ( ( val[0] >= val[1] ) || (val[0] <= 0.0) || (val[1] <= 0.0) ) {
                        throw Scheme::bad_scheme(Scheme::BadRange);
                    }
                    ranges.push_back(val);
                }
            }

            file.close();

            if ( ranges.empty() ) {
                throw Scheme::bad_scheme(Scheme::BadRange);
            }

        } catch (...) {
            file.close();
            throw;
        }
    }

    return ranges;
}



        /*  Scheme class public methods  */


// load and parse input scheme description (in XML form)
void Scheme::load_from_file(const string &filename, bool QE_flag)
{
    FILE *scheme_file;

    string sclass,stype,sshape,scomment;
    string bpr,btype,bshape,brange_file,brange_distr;
    string sname,resfileroot,logfile;

    Beam::BeamType beam_type;
    Beam::BeamShape beam_shape;
    Beam::BeamProfile beam_prof;
    Beam::BeamRangeDistr beam_range_distr;

    Surface::SurfaceClass cl;
    Surface::SurfaceType tp;
    Surface::SurfaceShape sh;

    Grating::GratingType gtype;

    Surface *surface;

    vector<real_t> spars,ssize,scen,sdist,sangle,sorder,sconst;
    vector<real_t> bpars, bcenter, brange;

//#if defined(WIN32) && defined(_MSC_VER)
//    uint8_t number_of_grating = 0; // number of grating in scheme
//#endif

//#if defined(__GNUC__) && !defined(USING_MACOSX)
//    u_int8_t number_of_grating = 0; // number of grating in scheme
//#endif

    uint8_t number_of_grating = 0; // number of grating in scheme

    // delete data if exists
    if ( !scheme.empty() ) destroy();
//    if ( !scheme.empty() ) scheme.clear();
//    wavelenth_ranges.clear();

    SchemeFileName = filename;

    scheme_file = fopen(filename.c_str(),"r");
    if ( scheme_file == NULL ) {
        throw bad_scheme(Scheme::FileOpenFailure);
    }

#ifdef USING_LINUX
    scheme_tree = mxmlLoadFile(NULL, scheme_file, MXML_OPAQUE_CALLBACK);
#endif

#ifdef USING_MSVC // On Windows platform mxmlLoadFile does not work!!!
    char *xml_str = NULL;
    struct stat filestatus;

    int st = stat(filename.c_str(),&filestatus); // ask for file size in bytes
    if ( st ) {
        throw bad_scheme(Scheme::FileOpenFailure);
    }

    try {
        xml_str = new char[filestatus.st_size+1];
        size_t n = fread(xml_str,sizeof(char),filestatus.st_size,scheme_file);
//        cout << "NNN===" << n << "; STAT: " << filestatus.st_size << endl;
//        if ( n != filestatus.st_size ) throw 1;
        xml_str[filestatus.st_size] = '\0';
//        xml_str[n] = '\0';

        scheme_tree = mxmlLoadString(NULL,xml_str,MXML_OPAQUE_CALLBACK);
    } catch (bad_alloc &ex) {
        fclose(scheme_file);
        throw bad_scheme(Scheme::MemAllocFailure);
    } catch (int ex) {
        delete[] xml_str;
        fclose(scheme_file);
        throw bad_scheme(Scheme::FileOpenFailure);
    }
    delete[] xml_str;
#endif

    fclose(scheme_file);

    if ( scheme_tree == NULL ) {
        cerr << "XML PARSING FAILURE!\n";
        throw bad_scheme(Scheme::ParseFailure);
    } else {
//        cerr << "XML PARSING OK!\n";
    }

    if ( QE_flag ) { // just load XML-tree for QE computation
        return;
    }

    // START SCHEME PARSING ...

    // read scheme general description
    desc_index = mxmlIndexNew(scheme_tree, "general", NULL);
    if ( desc_index == NULL ) {
        cerr << "BAD GENERAL DESCRIPTION INDEX\n";
        mxmlDelete(scheme_tree);
        throw bad_scheme(Scheme::ParseFailure);
    }

    // read the only first occurence of general scheme description
    current_node = mxmlIndexFind(desc_index, "general", NULL);
    if ( current_node == NULL ) {
        mxmlIndexDelete(desc_index);
        mxmlDelete(scheme_tree);
        throw bad_scheme(Scheme::NoSchemeDesc);
    }

    // mandatory attributes
    try {
        sname = get_string_attr("name",false);
        resfileroot = get_string_attr("result_file",false);
    } catch ( Scheme::bad_scheme &ex) { // rethrow exceptions
        mxmlIndexDelete(desc_index);
        mxmlDelete(scheme_tree);
        throw;
    }
    mxmlIndexDelete(desc_index);

    SchemeName = sname;
    ResultFileRoot = resfileroot;

    // optional log file attribute
    try {
        logfile = get_string_attr("log_file",false);
        LogFileName = logfile;
    } catch ( Scheme::bad_scheme &ex) { // just no log file
        LogFileName = "";
    }

    // read surface attributes and create list of surfaces (scheme) in memory

    surface_index = mxmlIndexNew(scheme_tree, "surface", NULL);
    if ( surface_index == NULL ) {
        cerr << "BAD INDEX\n";
        mxmlDelete(scheme_tree);
        throw bad_scheme(Scheme::ParseFailure);
    }
    mxmlIndexReset(surface_index);

    size_t i_surface = 0;

    while ( (current_node = mxmlIndexFind(surface_index, "surface", NULL)) != NULL ) {
//    while ( (current_node = mxmlIndexFind(surface_index, "surface", NULL)) != 0 ) {
//        cout << "SURF: " << i_surface << "; CN: " << current_node << endl;
        try {
            // first, read mandatory surface attributes

            // surface class

            sclass = get_string_attr("class");
//            cout << "CLASS: " << sclass << endl;
            if ( sclass != "AUX") { // special surface. only distance and angles are needed!
                // surface type
                stype = get_string_attr("type");

                // surface shape

                sshape = get_string_attr("shape");

                // surface parameters
                spars = get_numeric_attr("params");
//                if ( sclass != "PLANE" ) { // for plane one does not need parameters. it is just Z=0 plane
//                    spars = get_numeric_attr("params");
//                }

                // surface size and center

                ssize = get_numeric_attr("size");
                scen = get_numeric_attr("center");
            } else {
                stype = "AUX";
                sshape = "CIRC";
            }

            // distance, rotation angles and geometry

            sdist = get_numeric_attr("distance");
            sangle = get_numeric_attr("angles");

            // check mandatory attributes

            if ( sclass == "PLANE" ) {
                cl = Surface::Plane;
            } else if (sclass == "CONIC" ) {
                cl = Surface::Conic;
            } else if (sclass == "TORIC" ) {
                cl = Surface::Toric;
            } else if (sclass == "AUX") {
                cl = Surface::AuxClass;
            }
            else {
                throw bad_scheme(Scheme::BadAttrValue);
            }

            if ( stype == "MIRROR" ) {
                tp = Surface::Mirror;
            } else if ( stype == "LENS" ) {
                tp = Surface::Lens;
            } else if ( stype == "GRATING" ) {
                tp = Surface::Grating;
            } else if ( stype == "STOP" ) {
                tp = Surface::Stop;
            } else if ( stype == "SHADE" ) {
                tp = Surface::Shade;
            } else if ( stype == "DETECTOR" ) {
                tp = Surface::Detector;
            } else if ( stype == "AUX" ) {
                tp = Surface::AuxType;
            } else {
                throw bad_scheme(Scheme::BadAttrValue);
            }

            if ( !sshape.compare(0,4,"CIRC") ) {
                sh = Surface::Circle;
            } else if ( !sshape.compare(0,4,"RECT") ) {
                sh = Surface::Rectangle;
            } else {
                throw bad_scheme(Scheme::BadAttrValue);
            }

            // read optional comment
            try {
                scomment = get_string_attr("comment",false);
            } catch (Scheme::bad_scheme &ex) { // no comment attribute
                scomment = "";
            }

            // read specific surface attributes

            string sgtype;
            string s_n1_file;
            string s_n2_file;
            vector<real_t> n1;
            vector<real_t> n2;
            vector<long> sord;

            switch ( tp ) { // type of surface
                case Surface::Grating: { // reflective or transparent; for transparent: files/values with/of refractive indices as function of wavelength;
                                         // grating constant, working order or order range
                    try {
                        sgtype = get_string_attr("grating_type");
                        if ( sgtype == "REFLECTIVE") {
                            gtype = Grating::Reflective;
                        } else if ( sgtype == "TRANSPARENT" ) {
                            gtype = Grating::Transparent;
                        } else {
                            throw Grating::BadType;
                        }
                    } catch (Scheme::bad_scheme &ex) { // there is no "grating_type" attribute, interpretate it as default Reflective grating
                        gtype = Grating::Reflective;
                    }

                    sorder = get_numeric_attr("order");
                    // convert order to integer
                    long so;
                    for (size_t i = 0; i < sorder.size(); ++i) {
                        so = static_cast<long>(sorder[i]);
                        sord.push_back(so);
                    }

                    sconst = get_numeric_attr("const");

                    s_n1_file = "";
                    s_n2_file = "";
                    try {
                        n1 = get_numeric_attr("n1",true);
                    } catch (Scheme::bad_scheme &ex) { // it is a string, interpretate it as a filename
                        s_n1_file = get_string_attr("n1",false);
                    }
                    try {
                        n2 = get_numeric_attr("n2",true);
                    } catch (Scheme::bad_scheme &ex) { // it is a string, interpretating it as a filename
                        s_n2_file = get_string_attr("n2",false);
                    }

                    break;
                }
                case Surface::Lens: { // one needs files/values with/of refractive indices as function of wavelength
                    s_n1_file = "";
                    s_n2_file = "";
                    try {
                        n1 = get_numeric_attr("n1",true);
                    } catch (Scheme::bad_scheme &ex) { // it is a string, interpretate it as a filename
                        s_n1_file = get_string_attr("n1",false);
                    }
                    try {
                        n2 = get_numeric_attr("n2",true);
                    } catch (Scheme::bad_scheme &ex) { // it is a string, interpretating it as a filename
                        s_n2_file = get_string_attr("n2",false);
                    }
                    break;
                }
                case Surface::Mirror: { // no specific attributes
                    break;
                }
                case Surface::Stop: { // no specific attributes
                    break;
                }
                case Surface::Detector: { // no specific attributes
                    break;
                }
                case Surface::AuxType: { // no specific attributes
                    break;
                }
            }

            // create surface class instances
            switch ( tp ) {
                case Surface::Mirror: {
                    surface = new Mirror(cl,sh);
                    break;
                }
                case Surface::Lens: {
                    if ( s_n1_file.empty() ) {
                        if ( s_n2_file.empty() ) {
                            surface = new Lens(cl,sh,n1[0],n2[0]); // both ref. index are scalars
                        } else {
                            surface = new Lens(cl,sh,n1[0],s_n2_file.c_str());
                        }
                    } else {
                        if ( s_n2_file.empty() ) {
                            surface = new Lens(cl,sh,s_n1_file.c_str(),n2[0]);
                        } else {
                            surface = new Lens(cl,sh,s_n1_file.c_str(),s_n2_file.c_str());  // both ref. index are in files
                        }
                    }
                    break;
                }
                case Surface::Grating: {
                    if ( s_n1_file.empty() ) {
                        if ( s_n2_file.empty() ) {
                            surface = new Grating(cl,sh,gtype,sord,sconst[0],n1[0],n2[0]); // both ref. index are scalars
                        } else {
                            surface = new Grating(cl,sh,gtype,sord,sconst[0],n1[0],s_n2_file.c_str());
                        }
                    } else {
                        if ( s_n2_file.empty() ) {
                            surface = new Grating(cl,sh,gtype,sord,sconst[0],s_n1_file.c_str(),n2[0]);
                        } else {
                            surface = new Grating(cl,sh,gtype,sord,sconst[0],s_n1_file.c_str(),s_n2_file.c_str());  // both ref. index are in files
                        }
                    }
                    ++number_of_grating;
                    break;
                }
                case Surface::Stop: {
                    surface = new Stop(cl,sh);
                    break;
                }
                case Surface::Shade: {
                    surface = new Shade(cl,sh);
                    break;
                }
                case Surface::Detector: { // just generic class
                    surface = new Surface(cl,Surface::Detector,sh);
                    break;
                }
                case Surface::AuxType: {
                    surface = new Aux();
                    break;
                }
                default: break;
            }


            // define parameters
            surface->SetDistance(sdist);
            surface->SetAngles(sangle);
            surface->SetComment(scomment);

            if ( cl != Surface::AuxClass ) {
                surface->SetParams(spars);
                surface->SetGeometry(ssize,scen);
//                cout << "\nI_SURF: " << i_surface << "; CEN: " << scen[0] << "," << scen[1] << endl;
            }

            // add surface to scheme
            scheme.push_back(surface);

        } catch (bad_alloc &ex) {
            mxmlIndexDelete(surface_index);
            mxmlDelete(scheme_tree);
            throw bad_scheme(Scheme::MemAllocFailure);
        } catch (...) {
            mxmlIndexDelete(surface_index);
            mxmlDelete(scheme_tree);
            throw;
        }
        ++i_surface;
    }

    mxmlIndexDelete(surface_index);

    if ( i_surface == 0 ) {
        mxmlIndexDelete(surface_index);
        mxmlDelete(scheme_tree);

        throw bad_scheme(Scheme::NoSurfaceDesc);
//        cerr << "There were no surfaces in the input scheme desciption!!!\n";
    } else {
//        cerr << i_surface+1 << " surfaces are found.\n";
    }


    // read beam attributes and create beam in memory

    beam_index = mxmlIndexNew(scheme_tree, "beam", NULL);
    if ( beam_index == NULL ) {
        cerr << "BAD BEAM INDEX\n";
        mxmlDelete(scheme_tree);
        throw bad_scheme(Scheme::ParseFailure);
    }
    mxmlIndexReset(beam_index);

    // read the only first occurence of beam description
    current_node = mxmlIndexFind(beam_index, "beam", NULL);
    if ( current_node == NULL ) {
        mxmlIndexDelete(beam_index);
        mxmlDelete(scheme_tree);
        throw bad_scheme(Scheme::NoBeamDesc);
    }

    // beam mandatory attributes
    try {
        btype = get_string_attr("type");
        bshape = get_string_attr("shape");
        bpr = get_string_attr("profile");
        brange_distr = get_string_attr("range_distr");
        bpars = get_numeric_attr("params");
        bcenter = get_numeric_attr("center");
        wavelenth_ranges = get_range("range",number_of_grating);
    } catch ( Scheme::bad_scheme &ex) { // rethrow exceptions
        mxmlIndexDelete(beam_index);
        mxmlDelete(scheme_tree);
        throw;
    }

    mxmlIndexDelete(beam_index);

    // beam single range or file of ranges
    //
    // range must be a vector (or vectors if it is in file) of numbers and consists of:
    // w1 w2 [ord1 ord2 ...], where w1 and w2 are mandatory elements and are
    //                        minimal and maximal wavelength of each range;
    //                        optional ord1, ord2 etc are diffraction orders if
    //                        in the scheme gratings are given.
    //                        these (ord1, ord2 ...) elements are mandatory if
    //                        ranges are given in file and number gratings in
    //                        the scheme is greater than 0!!! If a single range is given
    //                        and there are gratings in the scheme then diffraction orders
    //                        are assumed to be given in the scheme file itself ("order" attribute of "Grating" surface).

//    bool isrange_file = true;
//    try {
//        brange_file = get_string_attr("range_file",false);
//    } catch ( Scheme::bad_scheme &ex)  {
//        isrange_file = false;
//    }

//    try {
//        brange = get_numeric_attr("range");
//    } catch ( Scheme::bad_scheme &ex)  {
//        if ( !isrange_file) { // there are neither file or single range. it is the error!
//            mxmlIndexDelete(beam_index);
//            mxmlDelete(scheme_tree);
//            throw Scheme::bad_scheme(Scheme::NoRange);
//        }
//    }
//    mxmlIndexDelete(beam_index);

//    // check range validity
//    AsciiFile range_file;
//    try {
//        if ( isrange_file ) { // try to open and read range file if it is given
//            //        ifstream range_file;
////            real_t w1, w2;
////            vector<real_t> vv(2,0.0);
//            vector<real_t> vv;

//            range_file.open(brange_file.c_str());

//            if ( range_file.fail() ) {
//                throw Scheme::bad_scheme(Scheme::NoRangeFile);
//            }

//            AsciiFile::AsciiFileFlag line_flag;
////            while ( (line_flag = range_file.ReadLine(2,&w1,&w2)) != AsciiFile::Eof ) {
//            while ( (line_flag = range_file.ReadLine(2+number_of_grating,vv)) != AsciiFile::Eof ) {
//                if ( line_flag == AsciiFile::InvalidData ) throw Scheme::bad_scheme(Scheme::BadRange);
//                if ( line_flag == AsciiFile::DataString ) {
////                    if ( ( w1 >= w2 ) || (w1 <= 0.0) || (w2 <= 0.0) ) {
//                    if ( ( vv[0] >= vv[1] ) || (vv[0] <= 0.0) || (vv[1] <= 0.0) ) {
//                        throw Scheme::bad_scheme(Scheme::BadRange);
//                    }
//                    wavelenth_ranges.push_back(vv);
//                }
//            }

//            range_file.close();

//            if ( wavelenth_ranges.empty() ) {
//                throw Scheme::bad_scheme(Scheme::BadRange);
//            }
//        } else { // a single range is given
//            if ( brange_distr != "MONOCHROMATIC" ) {
//                if ( (brange.size() < 2) || (brange[0] >= brange[1]) || (brange[0] <= 0.0) || (brange[1] <= 0.0) ) {
//                    throw Scheme::bad_scheme(Scheme::BadRange);
//                }
//            } else {
//                if ( brange[0] <= 0.0 ) throw Scheme::bad_scheme(Scheme::BadRange);
//            }
//            wavelenth_ranges.push_back(brange);
//        }
//    } catch (...) {
//        mxmlDelete(scheme_tree);
//        range_file.close();
//        throw;
//    }


    // create Beam type instance, set its parameters but do not compute beam itself (coordinates, cosins, wavelength vectors)

    size_t nlambda;
    bool composite_beam;

    try {
        if ( btype == "PARALLEL" ) {
            beam_type = Beam::Parallel;
        } else if ( btype == "CONIC" ) {
            beam_type = Beam::Conic;
        } else if ( btype == "GAUSS" ) {
            beam_type = Beam::Gauss;
        } else throw Beam::bad_beam(Beam::InvalidParams);

        if ( !bshape.compare(0,4,"CIRC") ) {
            beam_shape = Beam::Circle;
        } else if ( !bshape.compare(0,4,"RECT") ) {
            beam_shape = Beam::Rectangle;
        } else throw Beam::bad_beam(Beam::InvalidParams);

        if ( bpr == "RANDOM" ) {
            beam_prof = Beam::Random;
        } else if ( bpr == "CONCENTRIC" ) {
            beam_prof = Beam::Concentric;
        } else throw Beam::bad_beam(Beam::InvalidParams);

        if ( brange_distr == "RANDOM" ) {
            beam_range_distr = Beam::RandomDistr;
        } else if ( brange_distr == "UNIFORM" ) {
            beam_range_distr = Beam::UniformDistr;
        } else if ( brange_distr == "MONOCHROMATIC" ) {
            beam_range_distr = Beam::Monochromatic;
        } else throw Beam::bad_beam(Beam::InvalidParams);

        // optional attribute
        if ( beam_range_distr != Beam::Monochromatic) {
            try {
                vector<real_t> nl = get_numeric_attr("nlambda",true);
                nlambda = static_cast<size_t>(nl[0]);
                composite_beam = true;
            } catch (Scheme::bad_scheme &ex) {
                composite_beam = false;
            }
        } else {
            nlambda = 1;
            composite_beam = false;
        }

        beam = new Beam(beam_type,beam_shape,beam_prof,beam_range_distr);

        beam->SetParams(bpars,composite_beam,nlambda);
        beam->SetCenter(bcenter);
        beam->SetRange(wavelenth_ranges.front());

//        beam->Create(bpars,bcenter,wavelenth_ranges.front());

    } catch (...) {
        mxmlDelete(scheme_tree);
        throw;
    }

    mxmlDelete(scheme_tree);
}



void Scheme::destroy()
{
    delete beam;
    beam = NULL;
    wavelenth_ranges.clear();
    SchemeName = "Unknown";
    ResultFileRoot = "";
    LogFileName = "";
    total_QE_curve.clear();

    while( !scheme.empty() ) {
        delete scheme.front();
        scheme.pop_front();
    }
}


void Scheme::run(list<string> &resultfile_list, ostream &log_file, const size_t start_surface)
{
    string log_str;
    list<string> info_str;
    stringstream ss;

    time_t rawtime,endtime;
    struct tm *timeinfo;
    double elapsed_time;


    if ( scheme.empty() ) {
        log_file << "The scheme was not loaded!\n";
        throw Scheme::EmptyScheme;
    }

    if ( start_surface > scheme.size() ) {
        throw Scheme::BadStartIndex;
    }


    ss << "\nRAY-TRAYCING PACKAGE:  version " << RT_PACKAGE_MAJOR_VERSION << "." <<
          RT_PACKAGE_MINOR_VERSION << "\n";
    log_str = ss.str();
    log_file << log_str;

    log_str = "\n  RAY-TRACING ENGINE INFO:\n";
    log_file << log_str;

    // get engine info
    RT_engine_error err = rt_engine_info(info_str);
    if ( err == ENGINE_ERROR_OK ) {
        list<string>::iterator it;
//        for (size_t i = 0; i < info_str.size(); ++i ) {
        for (it = info_str.begin(); it != info_str.end(); ++it ) {
//            log_file << "    " << info_str[i] << endl;
            log_file << "    " << *it << endl;
        }
    }

    log_str = "\n  The scheme was loaded from: " + SchemeFileName + "\n";
    log_file << log_str;
    log_str = "  Scheme description: " + SchemeName;
    log_file << log_str << endl;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    resultfile_list.clear();

    ss.str("");
    string date = asctime(timeinfo);
    ss << "\nStarting ray-tracing procedure (" << date.substr(0,date.length()-1) << ") ...\n";
    log_str = ss.str();
    log_file << log_str;



    // start ray-tracing procedure
    size_t i_surface;
    size_t i_range;

    size_t i_grating;
    long diff_order;

    list<Surface*>::iterator current_surface;
    list<vector<real_t> >::iterator current_range;

    vector<real_t> range;

    vector<real_t> beam_dim;

    try {
        // save input beam
        beam->Recreate();
        ss.str("");
        ss << ResultFileRoot << "_input_beam.dat";
        beam->Save(ss.str());

        ss.str("");
        ss << "\n  Number of wavelength ranges: " << wavelenth_ranges.size() << "\n";
        log_file << ss.str();

        ss.str("");
        ss << "\n  Total number of rays in beam: " << beam->GetTotalNumberOfRays() << "\n";
        log_file << ss.str();

        current_range = wavelenth_ranges.begin();
        for ( i_range = 0; i_range < wavelenth_ranges.size(); ++i_range, ++current_range) {
            ss.str("");
            range = (*current_range);
            ss << "\n  *** Current range: [" << range[0] << ","  << range[1] << "] ***\n";
            log_file << ss.str();

            // create beam
            beam->SetRange(range);
            beam->Recreate();

            i_grating = 0;

            // skip "start_surface" surfaces
            current_surface = scheme.begin();
            for ( i_surface = 0; i_surface < start_surface; ++i_surface ) ++current_surface;

            // start ray-tracing cycle
            for ( i_surface = start_surface; i_surface < scheme.size(); ++i_surface,++current_surface ) {
                ss.str("");
                ss << "\n  Current surface index: " << i_surface << " ( " <<
                      (*current_surface)->GetComment() << ")\n";
                log_file << ss.str();

                log_file << "    Coordinate system transformation ...  ";
                beam->Transform(**current_surface);
                log_file << "OK\n";

                beam->Normalize();

                beam_dim = beam->GetBeamCosins();
                ss.str("");
                ss << "      Z-cosins range: [" << beam_dim[0] << "," << beam_dim[1] << "]" << endl;
                log_file << ss.str();

                if ( (*current_surface)->GetType() != Surface::AuxType ) {
                    log_file << "    Computation of beam-to-surface intersection coordinates ...  ";
                    (*current_surface)->Intersection(*beam);
                    ss.str("");
                    ss << "OK  ( number of good rays: " << beam->GetNumberOfGoodRays() << ")\n";
                    log_file << ss.str();

                    beam_dim = beam->GetBeamDim();
                    ss.str("");
                    ss << "      Incident ray dimension: x = [" << beam_dim[0] << "," << beam_dim[1] << "], y = [" <<
                          beam_dim[2] << "," << beam_dim[3] << "]" << endl;
                    log_file << ss.str();

                    log_file << "    Apply surface geometrical constrains ...  ";
                    (*current_surface)->ApplyConstrains(*beam);
                    ss.str("");
                    ss << "OK  ( number of good rays: " << beam->GetNumberOfGoodRays() << ")\n";
                    log_file << ss.str();

                    beam_dim = beam->GetBeamDim();
                    ss.str("");
                    ss << "      Ray dimension on the surface: x = [" << beam_dim[0] << "," << beam_dim[1] << "], y = [" <<
                          beam_dim[2] << "," << beam_dim[3] << "]" << endl;
                    log_file << ss.str();


                    if ( (*current_surface)->GetType() == Surface::Grating ) { // set diffraction order
                        if ( wavelenth_ranges.size() > 1 ) { // multiple orders! orders must be in range file
                            diff_order = static_cast<long>(range[2+i_grating]);
                            static_cast<Grating*>(*current_surface)->SetOrder(diff_order);
                            ++i_grating;
                        }
                    }

                    log_file << "    Surface action ...  ";
                    (*current_surface)->Action(*beam);
                    ss.str("");
                    ss << "OK  ( number of good rays: " << beam->GetNumberOfGoodRays() << ")\n";
                    log_file << ss.str();

                    beam->Normalize();

                    beam_dim = beam->GetBeamCosins();
                    ss.str("");
                    ss << "      Z-cosins range: [" << beam_dim[0] << "," << beam_dim[1] << "]" << endl;
                    log_file << ss.str();

                } else {

                }

//                ss.str("");
//                ss << ResultFileRoot << "_" << i_surface;
//                beam->Save(ss.str());

            }
            // save result beam
            ss.str("");
            if ( wavelenth_ranges.size() > 1 ) {
//                ss << ResultFileRoot << "_range_(" << range[0] << ":" << range[1] << ").dat";
                ss << ResultFileRoot << "_range_(" << range[0] << "-" << range[1] << ").dat";
            } else {
                ss << ResultFileRoot << ".dat";
            }

            beam->Save(ss.str());
            resultfile_list.push_back(ss.str());
        }
    } catch (...) {
        ss.str("");
        ss << "ERROR!\n";
        log_str = ss.str();
        log_file << log_str;
        throw;
    }

    time(&endtime);
    elapsed_time = difftime(endtime,rawtime); // in second

    unsigned long hours = 0;
    unsigned int minutes = 0;
    real_t seconds = elapsed_time;

    if ( elapsed_time > 3600.0 ) {
        hours =  static_cast<unsigned long>(seconds/3600.0); // hours
        seconds -= hours*3600.0;
    }
    minutes = static_cast<unsigned int>(seconds/60.0);
    seconds -= minutes*60.0;

    ss.str("");
//    ss << "\nFinished! ( elapsed time: " << elapsed_time << " minutes)\n";
    ss << "\nFinished! ( elapsed time: ";
    if ( hours ) ss << hours << "hours, ";
    ss << minutes << " minutes, " << seconds << " secs )\n";
    log_str = ss.str();
    log_file << log_str;

//    beam->Save(ResultFileRoot);
}


size_t Scheme::GetNumberOfSurfaces() const
{
    return scheme.size();
}


//
//  Compute quantum efficiency curve
//
//  inputs:
//    Nlambda - number of wavelength in each wavelength range where QE curve will be computed
//    QE_filename - name of result file
//
//  NOTE: The result file is binary unformatted one. It consists of blocks and
//  the number of the blocks is equal to a number of wavelength ranges.
//  Each block consists of a sequence of bytes in following order:
//     sizeof(size_t) - number of wavelength and QE values in the current block,
//     Nlambda*size(real_t) - wavelengths
//     Nlambda*sizeof(real_t) - QE values
//
//
void Scheme::ComputeQE(size_t Nlambda, char *QE_filename)
{
    // XML-tree already in scheme_tree (see load_from_file function)

    mxml_index_t *index;
    string sclass, stype, qe_file;
    vector<real_t> qe_val, val;
    real_t alpha, gamma, blaze_angle, grating_const;
    vector<long> order(1,0);
    TabulatedFunction qe_tab;
    Surface *surface;
    size_t i_surface = 0;

    uint8_t number_of_grating = 0; // number of grating in scheme

    try {
        // first, read surfaces description
        index = mxmlIndexNew(scheme_tree, "surface", NULL);
        if ( index == NULL ) {
            throw bad_scheme(Scheme::ParseFailure);
        }

        while ( (current_node = mxmlIndexFind(index, "surface", NULL)) != NULL ) {
            try { // "class" attribute is not mandatory now!
                sclass = get_string_attr("class");
            }  catch (Scheme::bad_scheme &ex) {
                sclass = "";
            }

            if ( sclass != "AUX") { // special surface. it will be skiped.
                // surface type
                stype = get_string_attr("type"); // it is mandatory attribute
            } else continue; // just read next surface

            if ( (stype == "MIRROR") || (stype == "LENS") || (stype == "GRATING") || (stype == "DETECTOR") ) {
                if ( stype == "GRATING" ) { // read grating attributes
                    val = get_numeric_attr("alpha");
                    alpha = val[0];
                    val = get_numeric_attr("gamma");
                    gamma = val[0];
                    val = get_numeric_attr("blaze_angle");
                    blaze_angle = val[0];
                    val = get_numeric_attr("const");
                    grating_const = val[0];
                    val = get_numeric_attr("order");
                    order[0] = static_cast<long>(val[0]);
                    surface = new Grating(Surface::Plane,Surface::Rectangle,Grating::Reflective,order,grating_const,1.0,1.0);
                    static_cast<Grating*>(surface)->SetRuleParams(blaze_angle,alpha,gamma);
//                    cout << "\n" << blaze_angle << ", " << alpha << ", " << gamma << ", " << order[0] << endl;
                    ++number_of_grating;
                } else {
                    surface = new Surface(Surface::Plane,Surface::Mirror,Surface::Circle);
                }

                // read surface "QE_curve" attribute

                try {
                    qe_val = get_numeric_attr("qe_curve",true); // a constant value is given
                    surface->SetQE(qe_val[0]);
                } catch (Scheme::bad_scheme &ex) { // a filename is given
                    if ( ex.scheme_error() == Scheme::ParseFailure ) throw; // no "qe_curve" attribute or it is bad
                    qe_file = get_string_attr("qe_curve",false);
                    surface->SetQE(qe_file.c_str());
                }

//                cout << "\nSTYPE: " << stype << ", QE file: " << qe_file << endl;

                scheme.push_back(surface);
                ++i_surface;
            } else {
                continue; // just read next surface
            }
        }

        mxmlIndexDelete(index);
        index = NULL;

        if ( !i_surface ) throw bad_scheme(Scheme::NoSurfaceDesc);

        // read beam description

        index = mxmlIndexNew(scheme_tree, "beam", NULL);
        if ( index == NULL ) {
            throw bad_scheme(Scheme::ParseFailure);
        }

        // read the only first occurence of beam description
        current_node = mxmlIndexFind(index, "beam", NULL);
        if ( current_node == NULL ) {
            throw bad_scheme(Scheme::NoBeamDesc);
        }


//cout << "range ..." << endl;

        // get "range" attribute
        wavelenth_ranges = get_range("range",number_of_grating);

//cout << "Nlambda: " << Nlambda << "; Nranges: " << wavelenth_ranges.size() << endl;

        mxmlIndexDelete(index);
        index = NULL;
        mxmlDelete(scheme_tree);
        scheme_tree = NULL;
    } catch (...) {
        mxmlIndexDelete(index);
        mxmlDelete(scheme_tree);
        throw;
    }

    // compute QE

    list<Surface*>::iterator current_surface = scheme.begin();
    list<vector<real_t> >::iterator current_range = wavelenth_ranges.begin();
    vector<real_t> lambda(Nlambda,0.0);
    real_t d_lambda;
    Surface::SurfaceType type;
    size_t vec_len, counts;

    FILE *file;

    file = fopen(QE_filename,"wb");
    if ( file == NULL ) {
        throw bad_scheme(Scheme::FileSaveFailure);
    }

    vec_len = sizeof(real_t)*Nlambda;

    uint8_t i_grating;

    try {
        for ( current_range = wavelenth_ranges.begin(); current_range != wavelenth_ranges.end(); ++current_range ) {
            val = (*current_range);

            d_lambda = (val[1]-val[0])/(Nlambda-1);

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    #ifdef USING_LINUX
            for (size_t i = 0; i < Nlambda; ++i ) {
    #endif
    #ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
            for (long long i = 0; i < Nlambda; ++i ) {
    #endif
                lambda[i] = val[0] + d_lambda*i;
            }

            qe_val = vector<real_t>(Nlambda,1.0); // init QE

            i_grating = 0;

//            cout << "\nrange: " << val[0] << ":" << val[1] << endl;

            for ( current_surface = scheme.begin(); current_surface != scheme.end(); ++current_surface) {
                type = (*current_surface)->GetType();
//                cout << "QE val: " << qe_val[50] << endl;
                switch ( type ) {
                    case Surface::Grating: {
                        if ( wavelenth_ranges.size() > 1 ) { // order from range file
                            order[0] = static_cast<long>(val[i_grating+2]);
                            static_cast<Grating*>(*current_surface)->SetOrder(order[0]);
                            ++i_grating;
                        }
                        static_cast<Grating*>(*current_surface)->ApplyQE(lambda,qe_val);
                        break;
                    }
                    case Surface::Mirror: ;
                    case Surface::Lens: {
                        (*current_surface)->ApplyQE(lambda,qe_val);
                        break;
                    }
                    default: ;
                }
//                cout << "TYPE: " << type << "; max: " << *max_element(qe_val.begin(),qe_val.end()) << "; min: " << *min_element(qe_val.begin(),qe_val.end()) << endl;
            }

            // save block

            counts = fwrite(&Nlambda,sizeof(size_t),1,file);
            if ( counts != 1 ) {
                throw bad_scheme(Scheme::FileSaveFailure);
            }

            counts = fwrite((void*)lambda.data(),sizeof(real_t),Nlambda,file);
            if ( counts != Nlambda ) {
                throw bad_scheme(Scheme::FileSaveFailure);
            }

            counts = fwrite((void*)qe_val.data(),sizeof(real_t),Nlambda,file);
            if ( counts != Nlambda ) {
                throw bad_scheme(Scheme::FileSaveFailure);
            }
        }
    } catch (...) {
        fclose(file);
        throw;
    }

    fclose(file);
}
