#ifndef ENUM_H
#define ENUM_H
typedef enum boundary {
	DIRICHLET = 0,
	NEUMANN = 3
} boundary;

typedef enum dMode {
	ONE = 0,
	FX1 = 8,
	FX2 = 1,
	CX2 = 7,
	BX1 = 6,
	BX2 = 2,
	FXX1 = 5,
	CXX2 = 4,
	BXX1 = 3
} dMode;

typedef enum dDir {
	X = 1,
	Y = 2,
	Z = 3,
	XY = 9,
	YZ = 8,
	ZX = 7
} dDir;

typedef enum side {
	LEFT = 2,
	RIGHT = 3,
	FRONT = 4,
	BACK = 5,
	SOUTH = 6,
	NORTH = 7,
	INTERIOR = 9
} Side;

#endif
