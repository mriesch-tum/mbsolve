#ifndef TYPES_H
#define TYPES_H

#include <complex>
#include <stdexcept>
#include <string>

namespace mbsolve {
/* type switch single/double */

/* complex number? boost? */

//typedef std::complex<double> complex;

/* use standard container */

typedef std::complex<double> complex;

typedef double real;
//typedef float real;


    /* TODO: make Result 2D? */
    /* TODO: Result assignment operator???? */
    /* TODO: use library for matrices? */


class Result
{
private:
    std::string m_name;
    unsigned int m_cols;
    unsigned int m_rows;
    unsigned int m_count;
    real *m_values;

    Result(const Result& other) { }

    Result& operator=(const Result& other) { return *this; }

public:
    explicit Result(const std::string& name, unsigned int cols,
		    unsigned int rows) :
	m_name(name), m_cols(cols), m_rows(rows), m_count(cols * rows)
    {
	m_values = new real[m_count];
    }

    ~Result() {
	delete[] m_values;
    }

    const std::string& name() const { return m_name; }

    unsigned int count() const { return m_count; }

    unsigned int cols() const { return m_cols; }

    unsigned int rows() const { return m_rows; }

    real *data(unsigned int row = 0) { return &m_values[row * m_cols]; }

    real& at(unsigned int index) {
	if (index > m_count) {
	    throw std::out_of_range("Index out of bounds");
	}
	return m_values[index];
    }

};

}

#endif
