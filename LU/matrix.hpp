#ifndef HEADER_MATRIX_HPP_INCLUDED
#define HEADER_MATRIX_HPP_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>

template <typename T>
class Matrix : private std::vector<T>
{
public:
    Matrix (const Matrix& other);
    Matrix (unsigned height = 0, unsigned width = 0);
    Matrix (unsigned height, unsigned width, const T& element);
    Matrix (unsigned height, unsigned width, T&& element);
    T &at(unsigned y, unsigned x);
    void fill (const T& element);
    void setHeight (unsigned height);
    void setWidth (unsigned width);
    void setSize (unsigned height, unsigned width);
    unsigned getHeight () const;
    unsigned getWidth () const;

    // Get iteration abilities from vector
    using std::vector<T>::begin;
    using std::vector<T>::end;
    using std::vector<T>::rbegin;
    using std::vector<T>::rend;
    using std::vector<T>::cbegin;
    using std::vector<T>::cend;
    using std::vector<T>::crbegin;
    using std::vector<T>::crend;

    // Get element management from vector
    using std::vector<T>::clear;
    using std::vector<T>::swap;

protected:
    bool isInRange (unsigned x, unsigned y) const;
private:
    unsigned height_, width_;

    friend std::ostream & operator << (std::ostream & ostr, Matrix & M)
    {
        ostr << "Matrix height: " << M.getHeight() << ", width: " << M.getWidth() << std::endl;
        for (unsigned i = 0; i < M.getHeight(); i++) {
            ostr << "        ";
            for (unsigned j = 0; j < M.getWidth(); j++)
                ostr << std::setw(4) << M.at(i, j) << " ";
            ostr << std::endl;
        }
        return ostr;
    }
    
    Matrix< T > & operator = (const Matrix< T > Other) const
    {
        std::vector<T>::operator =(Other);
        width_ = Other.width_;
        height_ = Other.height_;
    }

    Matrix< T > & operator = (const Matrix< T > & Other) const
    {
        std::vector<T>::operator =(Other);
        width_ = Other.width_;
        height_ = Other.height_;
    }

    Matrix< T > & operator = (const Matrix< T > && Other) const
    {
        std::vector<T>::operator =(Other);
        width_ = Other.width_;
        height_ = Other.height_;
    }
};

template <typename T>
Matrix<T>::Matrix(const Matrix& other) :
    std::vector<T>(other), height_(other.height_), width_(other.width_)
{
}

template <typename T>
Matrix<T>::Matrix (unsigned height, unsigned width):
    std::vector<T>(height * width, T()), height_ (height), width_ (width)
{

}

template <typename T>
Matrix<T>::Matrix (unsigned height, unsigned width, const T &element):
    std::vector<T>(height * width, T()), height_ (height), width_ (width)
{
    
}

template <typename T>
Matrix<T>::Matrix (unsigned height, unsigned width, T&& element):
    std::vector<T>(height * width, element), height_ (height), width_ (width)
{

}

template <typename T>
T &Matrix<T>::at(unsigned y, unsigned x)
{
    return std::vector<T>::at(y * width_ + x);
}


template <typename T>
bool Matrix<T>::isInRange(unsigned x, unsigned y) const
{
    return ((x < height_) && (y < width_));
}

template <typename T>
void Matrix<T>::fill(const T &element)
{
    std::fill(std::vector<T>::begin(), std::vector<T>::end(), element);
}

template <typename T>
unsigned Matrix<T>::getHeight() const
{
    return height_;
}

template <typename T>
unsigned Matrix<T>::getWidth() const
{
    return width_;
}

template <typename T>
void Matrix<T>::setHeight (unsigned height)
{
    setSize(height, width_);
}

template <typename T>
void Matrix<T>::setWidth(unsigned width)
{
    setSize(height_, width);
}

template <typename T>
void Matrix<T>::setSize(unsigned height, unsigned width)
{
    if (height == height_ && width == width_)
        return;
    
    unsigned lWidth = width_;
    unsigned lHeight = height_;
    Matrix< T > buffer(*this);

    width_ = width;
    height_ = height;

    std::vector<T>::resize(height_ * width_);
    std::fill(std::vector<T>::begin(), std::vector<T>::end(), T());

    for (unsigned i = 0; i < std::min(lHeight, height_); i++)
        for (unsigned j = 0; j < std::min(lWidth, width_); j++)
            at(i, j) = buffer.at(i, j);
}



#endif // MATRIX_HPP


