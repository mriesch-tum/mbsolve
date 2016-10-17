#ifndef QUANTITY_H
#define QUANTITY_H

#include <string>
#include <Types.hpp>

namespace mbsolve {

class Quantity
{
private:
    std::string m_name;
    real m_value;
public:
    Quantity(const real& value = 0.0);

    virtual const real& operator()() const;

    Quantity operator*(const Quantity& rhs) const;

};

static Quantity HBAR = 1.05457266e-34;
static Quantity MU0 =  M_PI * 4e-7;
static Quantity EPS0 = 8.854187817e-12;
static Quantity E0 =  1.60217733e-19;

}

#endif
