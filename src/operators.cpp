#include "operators.hpp"

std::unique_ptr<Field> operator+(Field & f1, double value)
{
	auto newField = std::make_unique<Field>(f1);

	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k) + value);
			}

	return newField;

}

std::unique_ptr<Field> operator+(Field & f1, Field & field)
{
	auto newField = std::make_unique<Field>(f1);

	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k) + field.get(i,j,k));
			}

	return newField;

}

std::unique_ptr<Field> operator-(Field & f1, Field & field)
{
	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k) - field.get(i,j,k));
			}

	return newField;

}

std::unique_ptr<Field> operator-(Field & f1, double value)
{
	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k) - value);
			}

	return newField;

}

std::unique_ptr<Field> operator/(Field & f1, Field & field)
{
	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k)/field.get(i,j,k));
			}

	return newField;
}

std::unique_ptr<Field> operator/(Field & f1, double value)
{
	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k)/value);
			}

	return newField;
}

std::unique_ptr<Field> operator/(double value, Field & f1)
{
	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, value/f1.get(i,j,k));
			}

	return newField;
}

std::unique_ptr<Field> operator^(Field & f1, double n)
{

	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, std::pow(f1.get(i,j,k), n));
			}

	return newField;
}

std::unique_ptr<Field> operator*(Field & f1, double factor)
{

	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k) * factor);
			}

	return newField;

}

std::unique_ptr<Field> operator*(Field & f1, Field & field)
{

	auto newField = std::make_unique<Field>(f1);
	for(int i=0; i<f1.nx; i++)
		for(int j=0; j<f1.ny; j++)
			for(int k=0; k<f1.nz; k++)
			{
				newField->set(i,j,k, f1.get(i,j,k) * field.get(i,j,k));
			}

	return newField;

}
