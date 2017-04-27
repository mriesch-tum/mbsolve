#ifndef QUANTITY_H
#define QUANTITY_H

#include <string>
#include <types.hpp>

namespace mbsolve {

class Quantity
{
private:
    std::string m_name;
    real m_value;
public:
    Quantity(const real& value = 0.0);

    virtual const real& operator()() const;

    operator real() const {
	return m_value;
    }

    /* TODO: operator with real? */

    Quantity operator+(const Quantity& rhs) const;

    Quantity& operator+=(const Quantity& rhs);

    Quantity operator-(const Quantity& rhs) const;

    Quantity& operator-=(const Quantity& rhs);

    Quantity operator*(const Quantity& rhs) const;

    Quantity& operator*=(const Quantity& rhs);

    Quantity operator/(const Quantity& rhs) const;

    Quantity& operator/=(const Quantity& rhs);

    bool operator<(const Quantity& rhs) const;

    bool operator>(const Quantity& rhs) const;

};

static Quantity HBAR = 1.05457266e-34;
static Quantity MU0 =  M_PI * 4e-7;
static Quantity EPS0 = 8.854187817e-12;
static Quantity E0 =  1.60217733e-19;

}

#endif
