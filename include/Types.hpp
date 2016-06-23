#ifndef TYPES_H
#define TYPES_H

#include <complex>
#include <string>

namespace mbsolve {
/* type switch single/double */

/* complex number? boost? */

typedef std::complex<double> complex;
typedef double real;


class Result
{
private:
    std::string m_name;
    unsigned int m_size;
    real *m_values;

    Result(const Result& other) { }

    Result& operator=(const Result& other) { return *this; }

public:
    explicit Result(const std::string& name) : m_name(name) { }

    ~Result() {
	delete[] m_values;
    }

    const std::string& name() { return m_name; }

    unsigned int size() { return m_size; }

    real *values() { return m_values; }

};

}

#endif
