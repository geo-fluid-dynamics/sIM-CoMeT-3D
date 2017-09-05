#ifndef OPERATORS_H
#define OPERATORS_H

#include <memory>
#include "field.hpp"

std::unique_ptr<Field> operator+(Field & f1, double value);
std::unique_ptr<Field> operator+(Field & f1, Field & field);
std::unique_ptr<Field> operator*(Field & f1, Field & field);
std::unique_ptr<Field> operator*(Field & f1, double factor);
std::unique_ptr<Field> operator^(Field & f1, double n);
std::unique_ptr<Field> operator/(Field & f1, double value);
std::unique_ptr<Field> operator/(Field & f1, Field & field);
std::unique_ptr<Field> operator/(double value, Field & f1);
std::unique_ptr<Field> operator-(Field & f1, Field & field);
std::unique_ptr<Field> operator-(Field & f1, double value);

std::unique_ptr<Field> operator+(std::unique_ptr<Field> f1, double value);

#endif
