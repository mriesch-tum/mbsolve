#include <Quantity.hpp>

namespace mbsolve {

Quantity::Quantity(const real& value) : m_value(value)
{
}

const real&
Quantity::operator()() const
{
    return m_value;
}

Quantity
Quantity::operator*(const Quantity& rhs) const
{
    return Quantity(this->m_value * rhs());
}


}
