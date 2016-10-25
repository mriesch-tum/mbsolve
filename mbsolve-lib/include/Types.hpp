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
    unsigned int m_size;
    real *m_values;

    Result(const Result& other) { }

    Result& operator=(const Result& other) { return *this; }

public:
    explicit Result(const std::string& name, unsigned int size) :
	m_name(name), m_size(size) {
	m_values = new real[size];
    }

    ~Result() {
	delete[] m_values;
    }

    const std::string& name() { return m_name; }

    unsigned int size() { return m_size; }

    real *data() { return m_values; }

    real& at(unsigned int index) {
	if (index > m_size) {
	    throw std::out_of_range("Index out of bounds");
	}
	return m_values[index];
    }

};

}

#endif
