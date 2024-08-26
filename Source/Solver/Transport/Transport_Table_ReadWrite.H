#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "Transport.H"

using namespace amrex;

template <typename TableType>
void c_TransportSolver::Read_Table1D(int assert_size, TableType &Tab1D_data,
                                     std::string filename)
{
    amrex::Print() << "Reading Table1D. filename: " << filename << "\n";

    std::ifstream infile;
    infile.open(filename.c_str());

    if (infile.fail())
    {
        amrex::Abort("Failed to read file " + filename);
    }
    else
    {
        auto const &Tab1D = Tab1D_data.table();
        auto thi = Tab1D_data.hi();
        auto tlo = Tab1D_data.lo();

        int filesize = 0;
        std::string line;
        while (infile.peek() != EOF)
        {
            std::getline(infile, line);
            filesize++;
        }
        amrex::Print() << "filesize: " << filesize << "\n";

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            filesize - 1 == assert_size,
            "Assert size, " + std::to_string(assert_size) +
                ", is not equal to the filesize-1, " +
                std::to_string(filesize - 1) + " !");

        infile.seekg(0, std::ios_base::beg);

        std::getline(infile, line);
        amrex::Print() << "file header: " << line << "\n";

        amrex::Real position, value;
        for (int l = tlo[0]; l < thi[0]; ++l)
        {
            infile >> position >> value;
            Tab1D(l) = value;
            amrex::Print() << "position/value: " << position << "    "
                           << Tab1D(l) << "\n";
        }
        infile.close();
    }
}

template <typename TableType>
void c_TransportSolver::Read_Table2D(int assert_size, TableType &Tab2D_data,
                                     std::string filename)
{
    amrex::Print() << "Reading Table2D. filename: " << filename << "\n";

    std::ifstream infile;
    infile.open(filename.c_str());

    if (infile.fail())
    {
        amrex::Abort("Failed to read file " + filename);
    }
    else
    {
        auto const &Tab2D = Tab2D_data.table();
        auto thi = Tab2D_data.hi();
        auto tlo = Tab2D_data.lo();

        int filesize = 0;
        std::string line;
        while (infile.peek() != EOF)
        {
            std::getline(infile, line);
            filesize++;
        }
        amrex::Print() << "filesize: " << filesize << "\n";

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            filesize - 1 == assert_size,
            "Assert size, " + std::to_string(assert_size) +
                ", is not equal to the filesize-1, " +
                std::to_string(filesize - 1) + " !");

        infile.seekg(0, std::ios_base::beg);

        std::getline(infile, line);
        amrex::Print() << "file header: " << line << "\n";

        amrex::Real value;
        for (int j = tlo[1]; j < thi[1]; ++j)  // slow moving index.
        {
            for (int i = tlo[0]; i < thi[0]; ++i)
            {
                infile >> value;
                Tab2D(i, j) = value;
            }
        }
        infile.close();
    }
}

template <typename VectorType, typename TableType>
void c_TransportSolver::Write_Table1D(const amrex::Vector<VectorType> &Vec,
                                      const TableType &Arr_data,
                                      std::string filename, std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        // amrex::Print() << "\nRoot Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const &Arr = Arr_data.const_table();
        auto thi = Arr_data.hi();

        outfile << header << "\n";
        if (Vec.size() == thi[0])
        {
            for (int e = 0; e < thi[0]; ++e)
            {
                outfile << std::setprecision(15) << std::setw(20) << Vec[e]
                        << std::setw(20) << Arr(e) << "\n";
            }
        }
        else
        {
            outfile << "Mismatch in the size of Vec and Table1D_data!"
                    << "\n";
        }
        outfile.close();
    }
}

template <typename TableType>
void c_TransportSolver::Write_Table2D(const TableData<TableType, 2> &Tab_data,
                                      std::string filename, std::string header)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        // amrex::Print() << "\nRoot Writing " << filename << "\n";
        std::ofstream outfile;
        outfile.open(filename.c_str());

        auto const &Tab2D = Tab_data.const_table();
        auto thi = Tab_data.hi();
        auto tlo = Tab_data.lo();

        outfile << header << "\n";
        for (int j = tlo[1]; j < thi[1]; ++j)  // slow moving index.
        {
            for (int i = tlo[0]; i < thi[0]; ++i)
            {
                outfile << std::setprecision(15) << Tab2D(i, j) << "\n";
            }
        }
        outfile.close();
    }
}
