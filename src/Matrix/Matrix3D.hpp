#ifndef MATRIX3D_HPP
#define MATRIX3D_HPP

#include "Matrix.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>      // ifstream
#include <sstream>      // istringstream
#include <stdexcept>    // std::runtime_error
#include <algorithm>    // std::equal()

// #define BUFFER_SIZE 4096
const size_t BUFFER_SIZE = 4096 ;

/*!
 * The Matrix3D class is a specialisation of the Matrix
 * class to make work with 3D matrices more easily.
 *
 * A text file format is defined to store such matrices. The specifications are as
 * follows :
 * Absolutely NO empty lines are allowed!
 * The following lines should contain :
 *
 * 1st line       : a slice header, ',,0' indicates that a slice of the 3rd dimension
 *                  is beginning (this is a z slice).
 * 2nd - Nth line : the firt slice, as a 2D matrix. In the example below, it has
 *                  dimensions 3x4.
 * N+1th line     : a slice header, ',,1' indicates that the 2nd slice is beginning.
 * N+1th - ...    : the second slice
 * and so on...
 *
 * Example of a 3x4x2 3D matrix
 * ---- start ----
 * ,,0
 *  1  2  3  4
 *  5  6  7  8
 *  8  9 10 11
 *,,1
 * 12 13 14 15
 * 16 17 18 19
 * 20 21 22 23
 * ----- end -----
 *
 */
template<class T>
class Matrix3D : public Matrix<T>
{
    public:
        // constructors
        Matrix3D() = default ;

        /*!
         * \brief Constructs a matrix with the given dimensions,
         * filled with 0 values.
         * \param dim1 the first dimension.
         * \param dim2 the second dimension.
         * \param dim3 the third dimension.
         */
        Matrix3D(size_t dim1, size_t dim2, size_t dim3) ;

        /*!
         * \brief Constructs a matrix with the given dimensions and
         * initialize the values to the given value.
         * \param dim1 the first dimension.
         * \param dim2 the second dimension.
         * \param dim3 the third dimension.
         * \param value the value to initialize the matrix content
         * with.
         */
        Matrix3D(size_t dim1, size_t dim2, size_t dim3, T value) ;

        /*!
         * \brief Copy constructor
         * \param other the matrix to copy the content from.
         */
        Matrix3D(const Matrix3D& other) ;

        /*!
         * \brief Constructs a matrix from a text file.
         * \param file_address the address of the file containing the matrix.
         * \throw std::runtime_error if anything happen while reading the
         * file (format error, file not found, etc).
         */
        Matrix3D(const std::string& file_address) throw (std::runtime_error) ;

        // methods overloaded from Matrix
        using Matrix<T>::get ;
        using Matrix<T>::set ;

        // methods
        /*!
         * \brief Gets the element at the given coordinates.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \throw std::out_of_range exception if the coordinates
         * are out of range.
         * \return the element.
         */
        T get(size_t dim1, size_t dim2, size_t dim3) const throw (std::out_of_range) ;
        /*!
         * \brief Sets the element at the given coordinates
         * to the given value.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \param value the new value.
         * \throw std::out_of_range exception if the coordinates
         * are out of range.
         */
        void set(size_t dim1, size_t dim2, size_t dim3, T value) throw (std::out_of_range) ;

        /*!
         * \brief Produces a nice representation of the matrix on the given
         * stream.
         * \param stream the stream.
         * \param precision the rounding precision.
         * \param width the column width in number of characters.
         * \param sep the character separator.
         */
        void print(std::ostream& stream, size_t precision=4 ,size_t width=6, char sep=' ') const ;

        // operators
        /*!
         * \brief Returns a reference to the corrresponding
         * element. This method does not perform any check on
         * the coordinates.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \return a reference to this element.
         */
        T& operator() (size_t dim1, size_t dim2, size_t dim3) ;
        /*!
         * \brief Returns a constant reference to the corrresponding
         * element. This method does not perform any check on
         * the coordinates.
         * \param dim1 the first dimension coordinate.
         * \param dim2 the second dimension coordinate.
         * \param dim3 the third dimension coordinate.
         * \return a constant reference to this element.
         */
        const T& operator() (size_t dim1, size_t dim2, size_t dim3) const ;

    private:
        // methods
        /*!
         * \brief Checks whether a given string is a slice header
         * (such as ",,0"), as found in files storing Matrix3D.
         * \param str the string to check.
         * \return whether the string is a slice header.
         */
        bool is_header(const std::string& str) const ;

} ;

/*!
 * \brief Sends a representation of the matrix to the stream.
 * \param stream the stream of interest.
 * \param m the matrix of interest.
 * \return a reference to the stream.
 */
template<class T>
std::ostream& operator << (std::ostream& stream, const Matrix3D<T>& m)
{   m.print(stream) ;
    return stream ;
}

template<class T>
Matrix3D<T>::Matrix3D(size_t dim1, size_t dim2, size_t dim3)
    : Matrix3D<T>(dim1, dim2, dim3, 0)
{}

template<class T>
Matrix3D<T>::Matrix3D(size_t dim1, size_t dim2, size_t dim3, T value)
    : Matrix<T>({dim1, dim2, dim3}, value)
{}

template<class T>
Matrix3D<T>::Matrix3D(const Matrix3D &other)
    : Matrix<T>(other)
{}


template<class T>
Matrix3D<T>::Matrix3D(const std::string &file_address) throw (std::runtime_error)
{
    this->_dim       = {0,0,0} ;
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
    std::vector<T>      buffer_vec ;
    T                   buffer_T ;

    // read file
    size_t n_line      = 0, n_line_data = 0 ; // number of line and of data line read
    size_t row_len     = 0, col_len = 0 ;     // length of row and column in nber of values
    size_t row_len_cur = 0, col_len_cur = 0 ; // current number of values read in row and col

    while(getline(file, buffer_str))
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
        // check whether it is the beginning of a slice
        // 1st line in file should be
        if(this->is_header(buffer_str))
        {   // check that slice have a constant number of rows
            if(this->_dim[2] == 1)
            {   col_len = col_len_cur ;
                this->_dim[0] = row_len ;
                this->_dim[1]  = col_len ;
            }
            else if(col_len_cur != col_len)
            {   file.close() ;
                char msg[BUFFER_SIZE] ;
                sprintf(msg, "format error! slice have variable dimensions in %s", file_address.c_str()) ;
                throw std::runtime_error(msg) ;
            }
            this->_dim[2]++ ;
            col_len_cur = 0 ;
            n_line++ ;
            continue ;
        }
        // 1st line in file should be a header and entering
        // this block is forbidden
        if(n_line == 0)
        {   file.close() ;
            char msg[BUFFER_SIZE] ;
            sprintf(msg, "format error! first line is not a slice header in %s", file_address.c_str()) ;
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
        {   file.close() ;
            char msg[BUFFER_SIZE] ;
            sprintf(msg, "format error! could not read a line in %s (incompatible data types)", file_address.c_str()) ;
            throw std::runtime_error(msg) ;
        }

        // check that number of column is constant
        if(n_line_data == 0)
        {  row_len = row_len_cur ; }
        else if(row_len_cur != row_len)
        {   file.close() ;
            char msg[BUFFER_SIZE] ;
            sprintf(msg, "format error! slice have variable dimensions in %s", file_address.c_str()) ;
            throw std::runtime_error(msg) ;
        }

        // update matrix content
        for(auto i : buffer_vec)
        {   this->_data.push_back(i) ;
            this->_data_size++ ;
        }
        col_len_cur++ ;
        n_line_data++ ;
        n_line++ ;
    }
    // check dimensions of last slice
    if(col_len_cur != col_len)
    {   file.close() ;
        char msg[BUFFER_SIZE] ;
        sprintf(msg, "format error! slice have variable dimensions in %s", file_address.c_str()) ;
        throw std::runtime_error(msg) ;
    }

    file.close() ;
    this->compute_dim_product() ;
}


template<class T>
T Matrix3D<T>::get(size_t dim1, size_t dim2, size_t dim3) const throw(std::out_of_range)
{   try
    {  return this->get({dim1, dim2, dim3}) ; }
    catch(std::out_of_range& e)
    {   throw e ; }
}

template<class T>
void Matrix3D<T>::set(size_t dim1, size_t dim2, size_t dim3, T value) throw(std::out_of_range)
{   try
    {  return this->set({dim1, dim2, dim3}, value) ; }
    catch(std::out_of_range& e)
    {   throw e ; }
}


template<class T>
T& Matrix3D<T>::operator () (size_t dim1, size_t dim2, size_t dim3)
{   std::vector<size_t> coord = {dim2, dim1, dim3} ;
    return this->_data[this->convert_to_offset(coord)] ;
}


template<class T>
void Matrix3D<T>::print(std::ostream& stream, size_t precision, size_t width, char sep) const
{
    stream.setf(std::ios::left) ;
    std::vector<size_t> dim = this->get_dim() ;
    for(size_t z=0; z<dim[2]; z++)
    {   stream << ",," << z << std::endl ;
        for(size_t x=0; x<dim[0]; x++)
        {   for(size_t y=0; y<dim[1]; y++)
            {   stream << std::setprecision(precision) << std::setw(width) << (*this)(x,y,z) << sep ; }
            stream << std::endl ;
        }
    }
}


template<class T>
const T& Matrix3D<T>::operator () (size_t dim1, size_t dim2, size_t dim3) const
{   std::vector<size_t> coord = {dim2, dim1, dim3} ;
    return this->_data[this->convert_to_offset(coord)] ;
}


template<class T>
bool Matrix3D<T>::is_header(const std::string& str) const
{   if(str[0] == ',' and
            str[1] == ',' and
            str.find(',', 2) == std::string::npos)
   {   return true ; }
   return false ;
}

#endif // MATRIX3D_HPP
