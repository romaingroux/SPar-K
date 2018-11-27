#ifndef MATRIX4D_HPP
#define MATRIX4D_HPP

#include "Matrix.hpp"

#include <string>
#include <vector>
#include <stdexcept> // std::out_of_range
#include <iostream>
#include <iomanip>
#include <fstream>  // ifstream
#include <sstream>  // sstream


// #define BUFFER_SIZE 4096
const size_t BUFFER_SIZE = 4096 ;

/*!
 * The Matrix4D class is a specialisation of the Matrix
 * class to handle 4D matrices more easily.
 *
 * A text file format is defined to store such matrices.
 * The specifications are as follows :
 *
 * Absolutely NO empty lines are allowed!
 * The following lines should contain :
 *
 * 1st line       : a slice header ',,,0' indicating that a
 *                  slice of the 4th dimension is beginning.
 * 3nd - Nth line : the slice of the 4th dimension. It contains
 *                  slice in the 3rd dimension which are 2D
 *                  matrices separated by headers (',,0' and
 *                  ',,1' in the below example).
 * N+1th line     : ',,,1' indicating that the 2nd slice of the
 *                  4th dimension is beginning.
 * and so on...
 *
 * Example
 * ---- start ----
 * ,,,0
 * ,,0
 * 1      2      3
 * 4      5      6
 * ,,1
 * 7      8      9
 * 10     11     12
 * ,,,1
 * ,,0
 * 21     22     23
 * 24     25     26
 * ,,1
 * 27     28     29
 * 30     31     32
 * ----- end -----
 *
 */
template<class T>
class Matrix4D : public Matrix<T>
{
    public:
        // constructors
        Matrix4D() = default ;
        /*!
         * \brief Constructs a matrix with the given dimensions,
         * filled with 0 values.
         * \param dim1 the first dimension.
         * \param dim2 the second dimension.
         * \param dim3 the third dimension.
         * \param dim4 the fourth dimension.
         */
        Matrix4D(size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4) ;
        /*!
         * \brief Constructs a matrix with the given dimensions and
         * initialize the values to the given value.
         * \param dim1 the first dimension.
         * \param dim2 the second dimension.
         * \param dim3 the third dimension.
         * \param dim4 the fourth dimension.
         * \param value the value to initialize the matrix content
         * with.
         */
        Matrix4D(size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4,
                 T value) ;
        /*!
         * \brief Copy constructor
         * \param other the matrix to copy the content from.
         */
        Matrix4D(const Matrix4D& other) ;
        /*!
         * \brief Constructs a matrix from a text file.
         * \param file_address the address of the file containing the matrix.
         * \throw std::runtime_error if anything happen while reading the
         * file (format error, file not found, etc).
         */
        Matrix4D(const std::string& file_address) throw (std::runtime_error) ;

        // methods overloaded from Matrix
        using Matrix<T>::get ;
        using Matrix<T>::set ;

        // methods OK
        /*!
         * \brief Gets the element at the given coordinates.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \param dim4 the fourth dimension coordinate.
         * \throw std::out_of_range exception if the coordinates
         * are out of range.
         * \return the element.
         */
        T get(size_t dim1,
              size_t dim2,
              size_t dim3,
              size_t dim4) const throw (std::out_of_range) ;
        /*!
         * \brief Sets the element at the given coordinates
         * to the given value.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \param dim4 the fourth dimension coordinate.
         * \param value the new value.
         * \throw std::out_of_range exception if the coordinates
         * are out of range.
         */
        void set(size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4,
                 T value) throw (std::out_of_range) ;
        /*!
         * \brief Produces a nice representation of the matrix on the given
         * stream.
         * \param stream the stream.
         * \param precision the rounding precision.
         * \param width the column width in number of characters.
         * \param sep the character separator.
         */
        void print(std::ostream& stream,
                   size_t precision=4,
                   size_t width=6,
                   char sep=' ') const ;

        // operators OK
        /*!
         * \brief Returns a reference to the corrresponding
         * element. This method does not perform any check on
         * the coordinates.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \param dim4 the third dimension coordinate.
         * \return a reference to this element.
         */
        T& operator() (size_t dim1,
                       size_t dim2,
                       size_t dim3,
                       size_t dim4) ;
        /*!
         * \brief Returns a reference to the corrresponding
         * element. This method does not perform any check on
         * the coordinates.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \param dim4 the third dimension coordinate.
         * \return a reference to this element.
         */
        const T& operator() (size_t dim1,
                             size_t dim2,
                             size_t dim3,
                             size_t dim4) const ;

    private:
        // methods
        /*!
         * \brief Checks whether a given string is a 3D header
         * (such as ",,0"), as found in files storing Matrix4D.
         * \param str the string to check.
         * \return whether the string is such a slice header.
         */
        bool is_header_3d(const std::string& str) const ;
        /*!
         * \brief Checks whether a given string is a 4D header
         * (such as ",,,0"), as found in files storing Matrix4D.
         * \param str the string to check.
         * \return whether the string is such a slice header.
         */
        bool is_header_4d(const std::string& str) const ;

        /*!
         * \brief Routine to load 4D matrices from files.
         * This method reads from a std::ifstream object,
         * from the current pointer location until i) a 4D
         * header line is found (such as ',,,1') or ii) until
         * it cannot read anymore from the stream. All
         * data are pushed back into the data vector and
         * the dimensions of the data read are stored into
         * the dim vector (these data are actually a 3D
         * matrix). If the method returned because it
         * found another 4D header, it returns true, false
         * otherwise.
         * To read an entire 4D matrix from a file, simply
         * use this scheme : i) read the 1st 4D header
         * ii) call this function while it returns true.
         * \param file_name a reference to a string containing
         * the address of the file currently read (for exception
         * messages).
         * \param file a reference to the std::ifstream to read
         * from. Obviously, the stream state will be modified as
         * the method reads from it. However, it will never be
         * closed by the method.
         * \param data a reference to an empty vector where the
         * read data will be pushed back.
         * \param dim a reference to an empty vector where the
         * dimensions of the read data will be stored.
         * \return whether the last piece of data read from the
         * stream was a 4D header.
         */
        bool get_3d_slice(const std::string& file_name,
                          std::ifstream& file,
                          std::vector<T>& data,
                          std::vector<size_t>& dim) const
        throw (std::runtime_error) ;

} ;

template<class T>
std::ostream& operator << (std::ostream& stream, const Matrix4D<T>& m)
{   m.print(stream) ;
    return stream ;
}

template<class T>
Matrix4D<T>::Matrix4D(size_t dim1,
                      size_t dim2,
                      size_t dim3,
                      size_t dim4)
    : Matrix<T>({dim1, dim2, dim3, dim4}, 0)
{}

template<class T>
Matrix4D<T>::Matrix4D(size_t dim1,
                      size_t dim2,
                      size_t dim3,
                      size_t dim4,
                      T value)
    : Matrix<T>({dim1, dim2, dim3, dim4}, value)
{}

template<class T>
Matrix4D<T>::Matrix4D(const Matrix4D &other)
    : Matrix<T>(other)
{}

template<class T>
Matrix4D<T>::Matrix4D(const std::string &file_address) throw (std::runtime_error)
{   this->_dim       = {0,0,0,0} ;
    this->_data      = std::vector<T>() ;
    this->_dim_size  = this->_dim.size() ;
    this->_data_size = this->_data.size() ;
    this->_dim_prod  = std::vector<size_t>(this->_dim_size, 0) ;

    std::ifstream file(file_address, std::ifstream::in) ;
    if(file.fail())
    {   char msg[BUFFER_SIZE] ;
        sprintf(msg, "error! cannot open %s", file_address.c_str()) ;
        throw std::runtime_error(msg) ;
    }

    std::string         buffer_str ;
    std::vector<T>      buffer_t ;
    std::vector<size_t> dim ;

    getline(file, buffer_str) ;
    bool found_4d_header = this->is_header_4d(buffer_str) ;
    do
    {   // check stream status and read content
        if(buffer_str.size() == 0)
        {   file.close() ;
            char msg[BUFFER_SIZE] ;
            sprintf(msg, "error! while reading %s (empty line)", file_address.c_str()) ;
            throw std::runtime_error(msg) ;
        }
        if(file.fail())
        {   file.close() ;
            char msg[BUFFER_SIZE] ;
            sprintf(msg, "error! while reading %s", file_address.c_str()) ;
            throw std::runtime_error(msg) ;
        }
        // this is the beginning of a 3D slice -> get it using routine
        if(found_4d_header)
        {   try
            {   // get slice
                buffer_t.clear() ;
                dim.clear() ;
                found_4d_header = this->get_3d_slice(file_address, file, buffer_t, dim);
                // update data
                for(const auto& i : buffer_t)
                {   this->_data.push_back(i) ;
                    this->_data_size++ ;
                }
                // update dim only for the 1st slice (the 1st slice set the dimensions)
                if(this->_dim[3] == 0)
                {   this->_dim[0] = dim[0] ;
                    this->_dim[1] = dim[1] ;
                    this->_dim[2] = dim[2] ;
                }
                // check dimensions of the slice
                else
                {   if(dim[0] != this->_dim[0] or
                       dim[1] != this->_dim[1] or
                       dim[2] != this->_dim[2])
                    {   char msg[BUFFER_SIZE] ;
                        sprintf(msg, "format error! slice have variable dimensions in %s",
                                file_address.c_str()) ;
                        throw std::runtime_error(msg) ;
                    }
                }
                this->_dim[3]++ ;
            }
            catch(std::runtime_error& e)
            {   file.close() ;
                throw e ;
            }
        }
    } while(found_4d_header) ;

    file.close() ;
    this->compute_dim_product() ;
}

template<class T>
T Matrix4D<T>::get(size_t dim1,
                   size_t dim2,
                   size_t dim3,
                   size_t dim4) const throw (std::out_of_range)
{   try
    {   return this->get({dim1, dim2, dim3, dim4}) ; }
    catch(std::out_of_range& e)
    {   throw e ; }
}

template<class T>
void Matrix4D<T>::set(size_t dim1,
                      size_t dim2,
                      size_t dim3,
                      size_t dim4,
                      T value) throw (std::out_of_range)
{   try
    {   this->set({dim1, dim2, dim3, dim4}, value) ; }
    catch(std::out_of_range& e)
    {   throw e ; }
}

template<class T>
void Matrix4D<T>::print(std::ostream &stream,
                        size_t precision,
                        size_t width,
                        char sep) const
{   stream.setf(std::ios::left) ;
    std::vector<size_t> dim = this->get_dim() ;

    for(size_t dim4=0; dim4<dim[3]; dim4++)
    {   stream << ",,," << dim4 << std::endl ;
        for(size_t dim3=0; dim3<dim[2]; dim3++)
        {   stream << ",," << dim3 << std::endl ;
            for(size_t dim2=0; dim2<dim[0]; dim2++)
            {   for(size_t dim1=0; dim1<dim[1]; dim1++)
                {   stream << std::setprecision(precision)
                           << std::setw(width)
                           << (*this)(dim2,dim1,dim3,dim4)
                           << sep ;
                }
                stream << std::endl ;
            }
        }
    }
}


template<class T>
T& Matrix4D<T>::operator () (size_t dim1,
                             size_t dim2,
                             size_t dim3,
                             size_t dim4)
{   std::vector<size_t> coord = {dim2, dim1, dim3, dim4} ;
    return this->_data[this->convert_to_offset(coord)] ;
}

template<class T>
const T& Matrix4D<T>::operator () (size_t dim1,
                                   size_t dim2,
                                   size_t dim3,
                                   size_t dim4) const
{   std::vector<size_t> coord = {dim2, dim1, dim3, dim4} ;
    return this->_data[this->convert_to_offset(coord)] ;
}

template<class T>
bool Matrix4D<T>::is_header_3d(const std::string &str) const
{   if(str[0] == ',' and
       str[1] == ',' and
       str.find(',', 2) == std::string::npos)
    {   return true ; }
    return false ;
}

template<class T>
bool Matrix4D<T>::is_header_4d(const std::string &str) const
{   if(str[0] == ',' and
       str[1] == ',' and
       str[2] == ',' and
       str.find(',', 3) == std::string::npos)
    {   return true ; }
    return false ;
}

template<class T>
bool Matrix4D<T>::get_3d_slice(const std::string& file_name,
                               std::ifstream& file,
                               std::vector<T> &data,
                               std::vector<size_t> &dim) const throw (std::runtime_error)
{
    bool found_4d_header = false ; // the flag to return

    dim = {0,0,0} ;

    std::string         buffer_str ;
    std::vector<T>      buffer_vec ;
    T                   buffer_T ;

    size_t n_line      = 0, n_line_data = 0 ; // number of line and of data line read
    size_t row_len     = 0, col_len = 0 ;     // length of row and column in nber of values
    size_t row_len_cur = 0, col_len_cur = 0 ; // current number of values read in row and col

    while(getline(file, buffer_str))
    {   // std::cerr << "in  " << n_line << " : " << buffer_str << std::endl ;
        // check stream status and read content
        if(buffer_str.size() == 0)
        {   char msg[BUFFER_SIZE] ;
            sprintf(msg, "error! while reading %s (empty line)", file_name.c_str()) ;
            throw std::runtime_error(msg) ;
        }
        if(file.fail())
        {   char msg[BUFFER_SIZE] ;
            sprintf(msg, "error! while reading %s", file_name.c_str()) ;
            throw std::runtime_error(msg) ;
        }
        // check whether this is the beginning of a 4D slice header, if so
        // break ;
        if(this->is_header_4d(buffer_str))
        {   found_4d_header = true ;
            break ;
        }
        // check whether it is the beginning of a slice
        // 1st line in file should be
        if(this->is_header_3d(buffer_str))
        {   // check that slice have a constant number of rows
            if(dim[2] == 1)
            {   col_len = col_len_cur ;
                dim[0] = row_len ;
                dim[1]  = col_len ;
            }
            else if(col_len_cur != col_len)
            {   char msg[BUFFER_SIZE] ;
                sprintf(msg, "format error! slice have variable dimensions in %s",
                        file_name.c_str()) ;
                throw std::runtime_error(msg) ;
            }
            dim[2]++ ;
            col_len_cur = 0 ;
            n_line++ ;
            continue ;
        }
        // 1st line in file should be a header and entering
        // this block is forbidden
        if(n_line == 0)
        {   char msg[BUFFER_SIZE] ;
            sprintf(msg, "format error! first line is not a slice header in %s",
                    file_name.c_str()) ;
            throw std::runtime_error(msg) ;
        }

        // parse line
        row_len_cur = 0 ;
        buffer_vec.clear() ;
        std::istringstream buffer_ss(buffer_str) ;
        while(buffer_ss >> buffer_T)
        {   buffer_vec.push_back(buffer_T) ;
            row_len_cur++ ;
        }
        // check for an error which likely indicates that a value could not be
        // casted into a type T (mixed data types in the file)
        if(buffer_ss.fail() and not buffer_ss.eof())
        {   char msg[BUFFER_SIZE] ;
            sprintf(msg, "format error! could not read a line in %s (incompatible data types)",
                    file_name.c_str()) ;
            throw std::runtime_error(msg) ;
        }

        // check that number of column is constant
        if(n_line_data == 0)
        {  row_len = row_len_cur ; }
        else if(row_len_cur != row_len)
        {   char msg[BUFFER_SIZE] ;
            sprintf(msg, "format error! slice have variable dimensions in %s",
                    file_name.c_str()) ;
            throw std::runtime_error(msg) ;
        }

        // update matrix content
        for(auto i : buffer_vec)
        {   data.push_back(i) ; }
        col_len_cur++ ;
        n_line_data++ ;
        n_line++ ;
    }

    // check dimensions of last slice
    if(col_len_cur != col_len)
    {   char msg[BUFFER_SIZE] ;
        sprintf(msg, "format error! slice have variable dimensions in %s",
                file_name.c_str()) ;
        throw std::runtime_error(msg) ;
    }

    return found_4d_header ;
}


#endif // MATRIX4D_HPP
