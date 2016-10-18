#include <Quantity.hpp>

namespace mbsolve {

/* TODO: use boost::quantity? */

Quantity::Quantity(const real& value) : m_value(value)
{
}

const real&
Quantity::operator()() const
{
    return m_value;
}

Quantity
Quantity::operator+(const Quantity& rhs) const
{
    return Quantity(this->m_value + rhs());
}

Quantity&
Quantity::operator+=(const Quantity& rhs)
{
    m_value += rhs();
    return *this;
}

Quantity
Quantity::operator-(const Quantity& rhs) const
{
    return Quantity(this->m_value - rhs());
}

Quantity&
Quantity::operator-=(const Quantity& rhs)
{
    m_value -= rhs();
    return *this;
}

Quantity
Quantity::operator*(const Quantity& rhs) const
{
    return Quantity(this->m_value * rhs());
}

Quantity&
Quantity::operator*=(const Quantity& rhs)
{
    m_value *= rhs();
    return *this;
}

Quantity
Quantity::operator/(const Quantity& rhs) const
{
    return Quantity(this->m_value / rhs());
}

Quantity&
Quantity::operator/=(const Quantity& rhs)
{
    m_value /= rhs();
    return *this;
}

bool
Quantity::operator<(const Quantity& rhs) const
{
    return this->m_value < rhs();
}

bool
Quantity::operator>(const Quantity& rhs) const
{
    return this->m_value > rhs();
}

}
