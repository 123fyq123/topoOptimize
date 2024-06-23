#pragma once
#include<Eigen/Core>
#include<map>
#include<fstream>
#include<iomanip>
#include<vector>
#include<iostream>
const int PlotPrecision = 11;
class OutputParaview {

public:
	void static ParaviewMulti(const Eigen::MatrixXd& points, std::map<std::string, Eigen::MatrixXd>& data,
        const std::vector<int>& np, const std::string& fn)
	{
        const int n = points.rows();

        std::string mfn(fn);
        mfn.append(".vts");
        std::ofstream file(mfn.c_str());
        file << std::fixed; // no exponents
        file << std::setprecision(PlotPrecision);

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
        file << "<StructuredGrid WholeExtent=\"0 " << np[0] - 1 << " 0 " << np[1] - 1 << " 0 "
            << (np.size() > 2 ? np[2] - 1 : 0) << "\">\n";
        file << "<Piece Extent=\"0 " << np[0] - 1 << " 0 " << np[0] - 1 << " 0 "
            << (np.size() > 2 ? np[0] - 1 : 0) << "\">\n";

        file << "<PointData>\n";
        for (typename std::map<std::string, Eigen::MatrixXd >::iterator it = data.begin(); it != data.end(); it++)
        {
            file << "<DataArray type=\"Float32\" Name=\"" << it->first << "\" format=\"ascii\" NumberOfComponents=\"" << (it->second.rows() == 1 ? 1 : 3) << "\">\n";
            if (it->second.rows() == 1)
                for (size_t j = 0; j < it->second.cols(); ++j)
                    file << it->second(0,j) << " ";
            else
            {
                for (size_t j = 0; j < it->second.cols(); ++j)
                {
                    for (size_t i = 0; i != it->second.rows(); ++i)
                        file << it->second(i, j) << " ";
                    for (size_t i = it->second.rows(); i < 3; ++i)
                        file << "0 ";
                }
            }
            file << "</DataArray>\n";
        }
        file << "</PointData>\n";
        file << "<Points>\n";
        file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
        for (size_t j = 0; j < points.cols(); ++j)
        {
            for (size_t i = 0; i != n; ++i)
                file << points(i, j) << " ";
            for (size_t i = n; i < 3; ++i)
                file << "0 ";
        }
        file << "</DataArray>\n";
        file << "</Points>\n";
        file << "</Piece>\n";
        file << "</StructuredGrid>\n";
        file << "</VTKFile>\n";

        file.close();
        
	}

public:
	//static 
};


class ParaviewCollection {
        typedef std::string String;
    public:

        /// Constructor using a filename.
        ParaviewCollection(std::string const& fn)
            : mfn(fn), counter(0)
        {
            mfile << "<?xml version=\"1.0\"?>\n";
            mfile << "<VTKFile type=\"Collection\" version=\"0.1\">";
            mfile << "<Collection>\n";
        }

        /// Adds a part in the collection, with complete filename (including extension) \a fn
        void addPart(String const& fn)
        {
            //assert(fn.find_last_of(".") != String::npos, "File without extension");
            assert(fn.find_last_of(".") != String::npos);
            //assert(counter != -1, "Error: collection has been already saved.");
            assert(counter != -1);
            mfile << "<DataSet part=\"" << counter++ << "\" file=\"" << fn << "\"/>\n";
        }

        /// Adds a part in the collection, with filename \a fn with extension \a ext appended
        void addPart(String const& fn, String const& ext)
        {
            //assert(counter != -1, "Error: collection has been already saved.");
            assert(counter != -1);
            mfile << "<DataSet part=\"" << counter++ << "\" file=\"" << fn << ext << "\"/>\n";
        }

        /// Adds a part in the collection, with filename \a fni and extension \a ext appended
        void addPart(String const& fn, int i, String const& ext)
        {
            //assert(counter != -1, "Error: collection has been already saved.");
            assert(counter != -1);
            mfile << "<DataSet part=\"" << i << "\" file=\"" << fn << i << ext << "\"/>\n";
        }


        void addTimestep(String const& fn, int tstep, String const& ext)
        {
            mfile << "<DataSet timestep=\"" << tstep << "\" file=\"" << fn << ext << "\"/>\n";
        }

        void addTimestep(String const& fn, int part, int tstep, String const& ext)
        {
            mfile << "<DataSet part=\""
                << part << "\" timestep=\""
                << tstep << "\" file=\""
                << fn << part << ext << "\"/>\n";//<<"_"
        }

        /// Finalizes the collection by closing the XML tags, always call
        /// this function (once) when you finish adding files
        void save()
        {
            //assert(counter != -1, "Error: gsParaviewCollection::save() already called.");
            assert(counter != -1);
            mfile << "</Collection>\n";
            mfile << "</VTKFile>\n";

            mfn.append(".pvd");
            std::ofstream f(mfn.c_str());
            //assert(f.is_open(), "Error creating " << mfn);
            assert(f.is_open());
            f << mfile.rdbuf();
            f.close();
            mfile.str("");
            counter = -1;
        }

    private:
        /// Pointer to char stream
        std::stringstream mfile;

        /// File name
        String mfn;

        /// Counter for the number of parts (files) added in the collection
        int counter;

    private:
        // Construction without a filename is not allowed
        ParaviewCollection();
};