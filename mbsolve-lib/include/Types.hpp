#ifndef TYPES_H
#define TYPES_H

#include <complex>
#include <stdexcept>
#include <string>

namespace mbsolve {
/* type switch single/double */

/* complex number? boost? */

//typedef std::complex<double> complex;
typedef double real;


    /* TODO: make Result 2D? */
    /* TODO: Result assignment operator???? */
    /* TODO: use library for matrices? */


class Result
{
private:
    std::string m_name;
    unsigned int m_cols;
    unsigned int m_rows;
    unsigned int m_size;
    real *m_values;

    Result(const Result& other) { }

    Result& operator=(const Result& other) { return *this; }

public:
    explicit Result(const std::string& name, unsigned int cols,
		    unsigned int rows) :
	m_name(name), m_cols(cols), m_rows(rows),
	m_size(cols * rows)
    {
	m_values = new real[rows * cols];
    }

    ~Result() {
	delete[] m_values;
    }

    const std::string& name() const { return m_name; }

    unsigned int size() const { return m_size; }

    unsigned int cols() const { return m_cols; }

    unsigned int rows() const { return m_rows; }

    real *data(unsigned int row = 0) { return &m_values[row * m_cols]; }

    real& at(unsigned int index) {
	if (index > m_size) {
	    throw std::out_of_range("Index out of bounds");
	}
	return m_values[index];
    }

};

}

#endif
