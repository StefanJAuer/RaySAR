/*******************************************************************************
 * matrices.cpp
 *
 * This module contains code to manipulate 4x4 matrices.
 *
 * ---------------------------------------------------------------------------
 * Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
 * Copyright 1991-2013 Persistence of Vision Raytracer Pty. Ltd.
 *
 * POV-Ray is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * POV-Ray is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------------
 * POV-Ray is based on the popular DKB raytracer version 2.12.
 * DKBTrace was originally written by David K. Buck.
 * DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
 * ---------------------------------------------------------------------------
 * $File: //depot/public/povray/3.x/source/backend/math/matrices.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/math/matrices.h"
#include "backend/math/vector.h"
#include "base/pov_err.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
* Local preprocessor defines
******************************************************************************/

/*  determinate |a b|
 *              |c d|  
 */
#define Det(a, b, c, d)      ((a)*(d) - (b)*(c))

/*****************************************************************************
* Local typedefs
******************************************************************************/



/*****************************************************************************
* Local variables
******************************************************************************/



/*****************************************************************************
* Static functions
******************************************************************************/



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   Initialize the matrix to the following values:
*
*     0.0   0.0   0.0   0.0
*     0.0   0.0   0.0   0.0
*     0.0   0.0   0.0   0.0
*     0.0   0.0   0.0   0.0
*
* CHANGES
*
*   -
*
******************************************************************************/

void MZero (MATRIX result)
{
	register int i, j;

	for (i = 0 ; i < 4 ; i++)
	{
		for (j = 0 ; j < 4 ; j++)
		{
			result[i][j] = 0.0;
		}
	}
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   Initialize the matrix to the following values:
*
*     1.0   0.0   0.0   0.0
*     0.0   1.0   0.0   0.0
*     0.0   0.0   1.0   0.0
*     0.0   0.0   0.0   1.0
*
* CHANGES
*
*   -
*
******************************************************************************/

void MIdentity (MATRIX result)
{
	register int i, j;

	for (i = 0 ; i < 4 ; i++)
	{
		for (j = 0 ; j < 4 ; j++)
		{
			if (i == j)
			{
				result[i][j] = 1.0;
			}
			else
			{
				result[i][j] = 0.0;
			}
		}
	}
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void MTimesA (MATRIX result, const MATRIX matrix2)
{
	DBL t0, t1, t2, t3;

	t0 = result[0][0];
	t1 = result[0][1];
	t2 = result[0][2];
	t3 = result[0][3];
	result[0][0] = t0 * matrix2[0][0] + t1 * matrix2[1][0] + t2 * matrix2[2][0] + t3 * matrix2[3][0];
	result[0][1] = t0 * matrix2[0][1] + t1 * matrix2[1][1] + t2 * matrix2[2][1] + t3 * matrix2[3][1];
	result[0][2] = t0 * matrix2[0][2] + t1 * matrix2[1][2] + t2 * matrix2[2][2] + t3 * matrix2[3][2];
	result[0][3] = t0 * matrix2[0][3] + t1 * matrix2[1][3] + t2 * matrix2[2][3] + t3 * matrix2[3][3];

	t0 = result[1][0];
	t1 = result[1][1];
	t2 = result[1][2];
	t3 = result[1][3];
	result[1][0] = t0 * matrix2[0][0] + t1 * matrix2[1][0] + t2 * matrix2[2][0] + t3 * matrix2[3][0];
	result[1][1] = t0 * matrix2[0][1] + t1 * matrix2[1][1] + t2 * matrix2[2][1] + t3 * matrix2[3][1];
	result[1][2] = t0 * matrix2[0][2] + t1 * matrix2[1][2] + t2 * matrix2[2][2] + t3 * matrix2[3][2];
	result[1][3] = t0 * matrix2[0][3] + t1 * matrix2[1][3] + t2 * matrix2[2][3] + t3 * matrix2[3][3];

	t0 = result[2][0];
	t1 = result[2][1];
	t2 = result[2][2];
	t3 = result[2][3];
	result[2][0] = t0 * matrix2[0][0] + t1 * matrix2[1][0] + t2 * matrix2[2][0] + t3 * matrix2[3][0];
	result[2][1] = t0 * matrix2[0][1] + t1 * matrix2[1][1] + t2 * matrix2[2][1] + t3 * matrix2[3][1];
	result[2][2] = t0 * matrix2[0][2] + t1 * matrix2[1][2] + t2 * matrix2[2][2] + t3 * matrix2[3][2];
	result[2][3] = t0 * matrix2[0][3] + t1 * matrix2[1][3] + t2 * matrix2[2][3] + t3 * matrix2[3][3];

	t0 = result[3][0];
	t1 = result[3][1];
	t2 = result[3][2];
	t3 = result[3][3];
	result[3][0] = t0 * matrix2[0][0] + t1 * matrix2[1][0] + t2 * matrix2[2][0] + t3 * matrix2[3][0];
	result[3][1] = t0 * matrix2[0][1] + t1 * matrix2[1][1] + t2 * matrix2[2][1] + t3 * matrix2[3][1];
	result[3][2] = t0 * matrix2[0][2] + t1 * matrix2[1][2] + t2 * matrix2[2][2] + t3 * matrix2[3][2];
	result[3][3] = t0 * matrix2[0][3] + t1 * matrix2[1][3] + t2 * matrix2[2][3] + t3 * matrix2[3][3];
}

void MTimesB (const MATRIX matrix1, MATRIX result)
{
	DBL t0, t1, t2, t3;

	t0 = result[0][0];
	t1 = result[1][0];
	t2 = result[2][0];
	t3 = result[3][0];
	result[0][0] = matrix1[0][0] * t0 + matrix1[0][1] * t1 + matrix1[0][2] * t2 + matrix1[0][3] * t3;
	result[1][0] = matrix1[1][0] * t0 + matrix1[1][1] * t1 + matrix1[1][2] * t2 + matrix1[1][3] * t3;
	result[2][0] = matrix1[2][0] * t0 + matrix1[2][1] * t1 + matrix1[2][2] * t2 + matrix1[2][3] * t3;
	result[3][0] = matrix1[3][0] * t0 + matrix1[3][1] * t1 + matrix1[3][2] * t2 + matrix1[3][3] * t3;

	t0 = result[0][1];
	t1 = result[1][1];
	t2 = result[2][1];
	t3 = result[3][1];
	result[0][1] = matrix1[0][0] * t0 + matrix1[0][1] * t1 + matrix1[0][2] * t2 + matrix1[0][3] * t3;
	result[1][1] = matrix1[1][0] * t0 + matrix1[1][1] * t1 + matrix1[1][2] * t2 + matrix1[1][3] * t3;
	result[2][1] = matrix1[2][0] * t0 + matrix1[2][1] * t1 + matrix1[2][2] * t2 + matrix1[2][3] * t3;
	result[3][1] = matrix1[3][0] * t0 + matrix1[3][1] * t1 + matrix1[3][2] * t2 + matrix1[3][3] * t3;

	t0 = result[0][2];
	t1 = result[1][2];
	t2 = result[2][2];
	t3 = result[3][2];
	result[0][2] = matrix1[0][0] * t0 + matrix1[0][1] * t1 + matrix1[0][2] * t2 + matrix1[0][3] * t3;
	result[1][2] = matrix1[1][0] * t0 + matrix1[1][1] * t1 + matrix1[1][2] * t2 + matrix1[1][3] * t3;
	result[2][2] = matrix1[2][0] * t0 + matrix1[2][1] * t1 + matrix1[2][2] * t2 + matrix1[2][3] * t3;
	result[3][2] = matrix1[3][0] * t0 + matrix1[3][1] * t1 + matrix1[3][2] * t2 + matrix1[3][3] * t3;

	t0 = result[0][3];
	t1 = result[1][3];
	t2 = result[2][3];
	t3 = result[3][3];
	result[0][3] = matrix1[0][0] * t0 + matrix1[0][1] * t1 + matrix1[0][2] * t2 + matrix1[0][3] * t3;
	result[1][3] = matrix1[1][0] * t0 + matrix1[1][1] * t1 + matrix1[1][2] * t2 + matrix1[1][3] * t3;
	result[2][3] = matrix1[2][0] * t0 + matrix1[2][1] * t1 + matrix1[2][2] * t2 + matrix1[2][3] * t3;
	result[3][3] = matrix1[3][0] * t0 + matrix1[3][1] * t1 + matrix1[3][2] * t2 + matrix1[3][3] * t3;
}

void MTimesC (MATRIX result, const MATRIX matrix1, const MATRIX matrix2)
{
	result[0][0] = matrix1[0][0] * matrix2[0][0] + matrix1[0][1] * matrix2[1][0] + matrix1[0][2] * matrix2[2][0] + matrix1[0][3] * matrix2[3][0];
	result[0][1] = matrix1[0][0] * matrix2[0][1] + matrix1[0][1] * matrix2[1][1] + matrix1[0][2] * matrix2[2][1] + matrix1[0][3] * matrix2[3][1];
	result[0][2] = matrix1[0][0] * matrix2[0][2] + matrix1[0][1] * matrix2[1][2] + matrix1[0][2] * matrix2[2][2] + matrix1[0][3] * matrix2[3][2];
	result[0][3] = matrix1[0][0] * matrix2[0][3] + matrix1[0][1] * matrix2[1][3] + matrix1[0][2] * matrix2[2][3] + matrix1[0][3] * matrix2[3][3];

	result[1][0] = matrix1[1][0] * matrix2[0][0] + matrix1[1][1] * matrix2[1][0] + matrix1[1][2] * matrix2[2][0] + matrix1[1][3] * matrix2[3][0];
	result[1][1] = matrix1[1][0] * matrix2[0][1] + matrix1[1][1] * matrix2[1][1] + matrix1[1][2] * matrix2[2][1] + matrix1[1][3] * matrix2[3][1];
	result[1][2] = matrix1[1][0] * matrix2[0][2] + matrix1[1][1] * matrix2[1][2] + matrix1[1][2] * matrix2[2][2] + matrix1[1][3] * matrix2[3][2];
	result[1][3] = matrix1[1][0] * matrix2[0][3] + matrix1[1][1] * matrix2[1][3] + matrix1[1][2] * matrix2[2][3] + matrix1[1][3] * matrix2[3][3];

	result[2][0] = matrix1[2][0] * matrix2[0][0] + matrix1[2][1] * matrix2[1][0] + matrix1[2][2] * matrix2[2][0] + matrix1[2][3] * matrix2[3][0];
	result[2][1] = matrix1[2][0] * matrix2[0][1] + matrix1[2][1] * matrix2[1][1] + matrix1[2][2] * matrix2[2][1] + matrix1[2][3] * matrix2[3][1];
	result[2][2] = matrix1[2][0] * matrix2[0][2] + matrix1[2][1] * matrix2[1][2] + matrix1[2][2] * matrix2[2][2] + matrix1[2][3] * matrix2[3][2];
	result[2][3] = matrix1[2][0] * matrix2[0][3] + matrix1[2][1] * matrix2[1][3] + matrix1[2][2] * matrix2[2][3] + matrix1[2][3] * matrix2[3][3];

	result[3][0] = matrix1[3][0] * matrix2[0][0] + matrix1[3][1] * matrix2[1][0] + matrix1[3][2] * matrix2[2][0] + matrix1[3][3] * matrix2[3][0];
	result[3][1] = matrix1[3][0] * matrix2[0][1] + matrix1[3][1] * matrix2[1][1] + matrix1[3][2] * matrix2[2][1] + matrix1[3][3] * matrix2[3][1];
	result[3][2] = matrix1[3][0] * matrix2[0][2] + matrix1[3][1] * matrix2[1][2] + matrix1[3][2] * matrix2[2][2] + matrix1[3][3] * matrix2[3][2];
	result[3][3] = matrix1[3][0] * matrix2[0][3] + matrix1[3][1] * matrix2[1][3] + matrix1[3][2] * matrix2[2][3] + matrix1[3][3] * matrix2[3][3];
}

/*  AAC - These are not used, so they are commented out to save code space...

void MAdd (MATRIX result, MATRIX matrix1, MATRIX matrix2)
{
	register int i, j;

	for (i = 0 ; i < 4 ; i++)
		for (j = 0 ; j < 4 ; j++)
			result[i][j] = (*matrix1)[i][j] + (*matrix2)[i][j];
}

void MSub (MATRIX result, MATRIX matrix1, MATRIX matrix2)
{
	register int i, j;

	for (i = 0 ; i < 4 ; i++)
		for (j = 0 ; j < 4 ; j++)
			result[i][j] = matrix1[i][j] - matrix2[i][j];
}

void MScale (MATRIX result, MATRIX matrix1, DBL amount)
{
	register int i, j;

	for (i = 0 ; i < 4 ; i++)
		for (j = 0 ; j < 4 ; j++)
			result[i][j] = matrix1[i][j] * amount;
}
... up to here! */



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

#define SWAP(a,b) c = a; a = b; b = c

void MTranspose (MATRIX result)
{
	DBL c;

	SWAP(result[0][1], result[1][0]);
	SWAP(result[0][2], result[2][0]);
	SWAP(result[1][2], result[2][1]);
	SWAP(result[2][3], result[3][2]);
	SWAP(result[3][0], result[0][3]);
	SWAP(result[3][1], result[1][3]);
}

#undef SWAP

void MTranspose (MATRIX result, const MATRIX matrix1)
{
	register int i, j;

	for (i = 0 ; i < 4 ; i++)
	{
		for (j = 0 ; j < 4 ; j++)
		{
			result[i][j] = matrix1[j][i];
		}
	}
}

/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Modified to not calculate anser_array[3]. [DB]
*
******************************************************************************/

void MTransPoint (VECTOR result, const VECTOR vector, const TRANSFORM *transform)
{
	register int i;
	DBL answer_array[4];
	const MATRIX *matrix;

	matrix = &transform->matrix;

	for (i = 0 ; i < 3 ; i++)
	{
		answer_array[i] = vector[X] * (*matrix)[0][i] +
		                  vector[Y] * (*matrix)[1][i] +
		                  vector[Z] * (*matrix)[2][i] + (*matrix)[3][i];
	}

	result[X] = answer_array[0];
	result[Y] = answer_array[1];
	result[Z] = answer_array[2];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Modified to not calculate anser_array[3]. [DB]
*
******************************************************************************/

void MInvTransPoint (VECTOR result, const VECTOR vector, const TRANSFORM *transform)
{
	register int i;
	DBL answer_array[4];
	const MATRIX *matrix;

	matrix = &transform->inverse;

	for (i = 0 ; i < 3 ; i++)
	{
		answer_array[i] = vector[X] * (*matrix)[0][i] +
		                  vector[Y] * (*matrix)[1][i] +
		                  vector[Z] * (*matrix)[2][i] + (*matrix)[3][i];
	}

	result[X] = answer_array[0];
	result[Y] = answer_array[1];
	result[Z] = answer_array[2];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Modified to not calculate anser_array[3]. [DB]
*
******************************************************************************/

void MTransDirection (VECTOR result, const VECTOR vector, const TRANSFORM *transform)
{
	register int i;
	DBL answer_array[4];
	const MATRIX *matrix;

	matrix = &transform->matrix;

	for (i = 0 ; i < 3 ; i++)
	{
		answer_array[i] = vector[X] * (*matrix)[0][i] +
		                  vector[Y] * (*matrix)[1][i] +
		                  vector[Z] * (*matrix)[2][i];
	}

	result[X] = answer_array[0];
	result[Y] = answer_array[1];
	result[Z] = answer_array[2];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   Sep 1994 : Modified to not calculate anser_array[3]. [DB]
*
******************************************************************************/

void MInvTransDirection (VECTOR result, const VECTOR vector, const TRANSFORM*transform)
{
	register int i;
	DBL answer_array[4];
	const MATRIX *matrix;

	matrix = &transform->inverse;

	for (i = 0 ; i < 3 ; i++)
	{
		answer_array[i] = vector[X] * (*matrix)[0][i] +
		                  vector[Y] * (*matrix)[1][i] +
		                  vector[Z] * (*matrix)[2][i];
	}

	result[X] = answer_array[0];
	result[Y] = answer_array[1];
	result[Z] = answer_array[2];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void MTransNormal (VECTOR result, const VECTOR vector, const TRANSFORM*transform)
{
	register int i;
	DBL answer_array[3];
	const MATRIX *matrix;

	matrix = &transform->inverse;

	for (i = 0 ; i < 3 ; i++)
	{
		answer_array[i] = vector[X] * (*matrix)[i][0] +
		                  vector[Y] * (*matrix)[i][1] +
		                  vector[Z] * (*matrix)[i][2];
	}

	result[X] = answer_array[0];
	result[Y] = answer_array[1];
	result[Z] = answer_array[2];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void MInvTransNormal (VECTOR result, const VECTOR vector, const TRANSFORM*transform)
{
	register int i;
	DBL answer_array[3];
	const MATRIX *matrix;

	matrix = &transform->matrix;

	for (i = 0 ; i < 3 ; i++)
	{
		answer_array[i] = vector[X] * (*matrix)[i][0] +
		                  vector[Y] * (*matrix)[i][1] +
		                  vector[Z] * (*matrix)[i][2];
	}

	result[X] = answer_array[0];
	result[Y] = answer_array[1];
	result[Z] = answer_array[2];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void Compute_Scaling_Transform (TRANSFORM *result, const VECTOR vector)
{
	MIdentity (result->matrix);

	(result->matrix)[0][0] = vector[X];
	(result->matrix)[1][1] = vector[Y];
	(result->matrix)[2][2] = vector[Z];

	MIdentity (result->inverse);

	(result->inverse)[0][0] = 1.0 / vector[X];
	(result->inverse)[1][1] = 1.0 / vector[Y];
	(result->inverse)[2][2] = 1.0 / vector[Z];
}



/*****************************************************************************
*
* FUNCTION
*
*   Compute_Matrix_Transform
*
* INPUT
*
*   matrix - matrix from which to create transform
*
* OUTPUT
*
*   result - complete transform
*
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
*
* DESCRIPTION
*
*   Builds a complete transform from a matrix.
*
* CHANGES
*
*   June 1995 : Creation
*
******************************************************************************/

void Compute_Matrix_Transform (TRANSFORM *result, const MATRIX matrix)
{
	register int i;

	for (i = 0; i < 4; i++)
	{
		(result->matrix)[i][0] = matrix[i][0];
		(result->matrix)[i][1] = matrix[i][1];
		(result->matrix)[i][2] = matrix[i][2];
		(result->matrix)[i][3] = matrix[i][3];
	}

	MInvers(result->inverse, result->matrix);
}



/* AAC - This is not used, so it's commented out...

void Compute_Inversion_Transform (TRANSFORM *result)
{
	MIdentity (result->matrix);

	(result->matrix)[0][0] = -1.0;
	(result->matrix)[1][1] = -1.0;
	(result->matrix)[2][2] = -1.0;
	(result->matrix)[3][3] = -1.0;


	(result->inverse)[0][0] = -1.0;
	(result->inverse)[1][1] = -1.0;
	(result->inverse)[2][2] = -1.0;
	(result->inverse)[3][3] = -1.0;
}
... up to here! */



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void Compute_Translation_Transform (TRANSFORM *transform, const VECTOR vector)
{
	MIdentity (transform->matrix);

	(transform->matrix)[3][0] = vector[X];
	(transform->matrix)[3][1] = vector[Y];
	(transform->matrix)[3][2] = vector[Z];

	MIdentity (transform->inverse);

	(transform->inverse)[3][0] = -vector[X];
	(transform->inverse)[3][1] = -vector[Y];
	(transform->inverse)[3][2] = -vector[Z];
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void Compute_Rotation_Transform (TRANSFORM *transform, const VECTOR vector)
{
	register DBL cosx, cosy, cosz, sinx, siny, sinz;
	MATRIX Matrix;
	VECTOR Radian_Vector;

	VScale (Radian_Vector, vector, M_PI_180);

	MIdentity (transform->matrix);

	cosx = cos (Radian_Vector[X]);
	sinx = sin (Radian_Vector[X]);
	cosy = cos (Radian_Vector[Y]);
	siny = sin (Radian_Vector[Y]);
	cosz = cos (Radian_Vector[Z]);
	sinz = sin (Radian_Vector[Z]);

	(transform->matrix) [1][1] = cosx;
	(transform->matrix) [2][2] = cosx;
	(transform->matrix) [1][2] = sinx;
	(transform->matrix) [2][1] = 0.0 - sinx;

	MTranspose (transform->inverse, transform->matrix);

	MIdentity (Matrix);

	Matrix [0][0] = cosy;
	Matrix [2][2] = cosy;
	Matrix [0][2] = 0.0 - siny;
	Matrix [2][0] = siny;

	MTimesA (transform->matrix, Matrix);

	MTranspose (Matrix);

	MTimesB (Matrix, transform->inverse);

	MIdentity (Matrix);

	Matrix [0][0] = cosz;
	Matrix [1][1] = cosz;
	Matrix [0][1] = sinz;
	Matrix [1][0] = 0.0 - sinz;

	MTimesA (transform->matrix, Matrix);

	MTranspose (Matrix);

	MTimesB (Matrix, transform->inverse);
}



/* AAC - This is not used so it's commented out...

void Compute_Look_At_Transform (TRANSFORM *result, VECTOR Look_At, VECTOR Up, VECTOR Right)
{
	MIdentity (result->inverse);

	(result->matrix)[0][0] = Right[X];
	(result->matrix)[0][1] = Right[Y];
	(result->matrix)[0][2] = Right[Z];

	(result->matrix)[1][0] = Up[X];
	(result->matrix)[1][1] = Up[Y];
	(result->matrix)[1][2] = Up[Z];

	(result->matrix)[2][0] = Look_At[X];
	(result->matrix)[2][1] = Look_At[Y];
	(result->matrix)[2][2] = Look_At[Z];

	MIdentity (result->matrix);

	MTranspose (result->matrix, result->inverse);
}
... up to here! */



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void Compose_Transforms (TRANSFORM *Original_Transform, const TRANSFORM *Additional_Transform)
{
	MTimesA(Original_Transform->matrix, Additional_Transform->matrix);

	MTimesB(Additional_Transform->inverse, Original_Transform->inverse);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*
* OUTPUT
*
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   Rotation about an arbitrary axis - formula from:
*
*     "Computational Geometry for Design and Manufacture", Faux & Pratt
*
*   NOTE: The angles for this transform are specified in radians.
*
* CHANGES
*
*   -
*
******************************************************************************/

void Compute_Axis_Rotation_Transform (TRANSFORM *transform, const VECTOR AxisVect, DBL angle)
{
	DBL cosx, sinx;
	VECTOR V1;

	VNormalize(V1, AxisVect);

	MIdentity(transform->matrix);

	cosx = cos(angle);
	sinx = sin(angle);

	transform->matrix[0][0] = V1[X] * V1[X] + cosx * (1.0 - V1[X] * V1[X]);
	transform->matrix[0][1] = V1[X] * V1[Y] * (1.0 - cosx) + V1[Z] * sinx;
	transform->matrix[0][2] = V1[X] * V1[Z] * (1.0 - cosx) - V1[Y] * sinx;

	transform->matrix[1][0] = V1[X] * V1[Y] * (1.0 - cosx) - V1[Z] * sinx;
	transform->matrix[1][1] = V1[Y] * V1[Y] + cosx * (1.0 - V1[Y] * V1[Y]);
	transform->matrix[1][2] = V1[Y] * V1[Z] * (1.0 - cosx) + V1[X] * sinx;

	transform->matrix[2][0] = V1[X] * V1[Z] * (1.0 - cosx) + V1[Y] * sinx;
	transform->matrix[2][1] = V1[Y] * V1[Z] * (1.0 - cosx) - V1[X] * sinx;
	transform->matrix[2][2] = V1[Z] * V1[Z] + cosx * (1.0 - V1[Z] * V1[Z]);

	MTranspose(transform->inverse, transform->matrix);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   Given a point and a direction and a radius, find the transform
*   that brings these into a canonical coordinate system
*
* CHANGES
*
*   7/24/95 Eduard Schwan  - Changed "if" condition to use EPSILON, not equality
*  12/12/95 Steve Demlow   - Clipped abs(up[Z]) to 1 to avoid acos overflow
*
******************************************************************************/

void Compute_Coordinate_Transform(TRANSFORM *trans, const VECTOR origin, VECTOR up, DBL radius, DBL length)
{
	TRANSFORM trans2;
	VECTOR tmpv;

	Make_Vector(tmpv, radius, radius, length);

	Compute_Scaling_Transform(trans, tmpv);

	if (fabs(up[Z]) > 1.0 - EPSILON)
	{
		Make_Vector(tmpv, 1.0, 0.0, 0.0);
		up[Z] = up[Z] < 0.0 ? -1.0 : 1.0;
	}
	else
	{
		Make_Vector(tmpv, -up[Y], up[X], 0.0);
	}

	Compute_Axis_Rotation_Transform(&trans2, tmpv, acos(up[Z]));

	Compose_Transforms(trans, &trans2);

	Compute_Translation_Transform(&trans2, origin);

	Compose_Transforms(trans, &trans2);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

TRANSFORM *Create_Transform()
{
	TRANSFORM *New;

	New = reinterpret_cast<TRANSFORM *>(POV_MALLOC(sizeof (TRANSFORM), "transform"));

	MIdentity (New->matrix);
	MIdentity (New->inverse);

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

TRANSFORM *Copy_Transform (const TRANSFORM*Old)
{
	TRANSFORM *New;
	if (Old != NULL)
	{
		New  = Create_Transform ();
		*New = *Old;
	}
	else
	{
		New = NULL;
	}

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

void Destroy_Transform (TRANSFORM *Trans)
{
	if(Trans != NULL)
		POV_FREE(Trans);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

VECTOR *Create_Vector ()
{
	VECTOR *New;

	New = reinterpret_cast<VECTOR *>(POV_MALLOC(sizeof (VECTOR), "vector"));

	Make_Vector (*New, 0.0, 0.0, 0.0);

	return (New);
}



/*****************************************************************************
*
* FUNCTION
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*   
* AUTHOR
*
*   POV-Ray Team
*   
* DESCRIPTION
*
*   -
*
* CHANGES
*
*   -
*
******************************************************************************/

DBL *Create_Float ()
{
	DBL *New_Float;

	New_Float = reinterpret_cast<DBL *>(POV_MALLOC(sizeof (DBL), "float"));

	*New_Float = 0.0;

	return (New_Float);
}



/*****************************************************************************
*
* FUNCTION
*
*   MInvers3
*
* INPUT
*
* OUTPUT
*
* RETURNS
*
* AUTHOR
*
*   Alexander Enzmann
*
* DESCRIPTION
*
*   - Invert a 3x3 Matrix
*
* CHANGES
*
*   -
*
******************************************************************************/
int MInvers3(const VECTOR inM[3], VECTOR  outM[3])
{
	DBL det;

	outM[0][X] =   (inM[1][Y] * inM[2][Z] - inM[1][Z] * inM[2][Y]);
	outM[1][X] = - (inM[0][Y] * inM[2][Z] - inM[0][Z] * inM[2][Y]);
	outM[2][X] =   (inM[0][Y] * inM[1][Z] - inM[0][Z] * inM[1][Y]);

	outM[0][Y] = - (inM[1][X] * inM[2][Z] - inM[1][Z] * inM[2][X]);
	outM[1][Y] =   (inM[0][X] * inM[2][Z] - inM[0][Z] * inM[2][X]);
	outM[2][Y] = - (inM[0][X] * inM[1][Z] - inM[0][Z] * inM[1][X]);

	outM[0][Z] =   (inM[1][X] * inM[2][Y] - inM[1][Y] * inM[2][X]);
	outM[1][Z] = - (inM[0][X] * inM[2][Y] - inM[0][Y] * inM[2][X]);
	outM[2][Z] =   (inM[0][X] * inM[1][Y] - inM[0][Y] * inM[1][X]);

	det = inM[0][X] * inM[1][Y] * inM[2][Z] +
	      inM[0][Y] * inM[1][Z] * inM[2][X] +
	      inM[0][Z] * inM[1][X] * inM[2][Y] -
	      inM[0][Z] * inM[1][Y] * inM[2][X] -
	      inM[0][X] * inM[1][Z] * inM[2][Y] -
	      inM[0][Y] * inM[1][X] * inM[2][Z];

	if (fabs(det) < 1.0e-10)
	{
		return (0);
	}

	det = 1.0 / det;

	VScaleEq(outM[0], det);
	VScaleEq(outM[1], det);
	VScaleEq(outM[2], det);

	return (1);
}



/*****************************************************************************
*
* FUNCTION
*
*   MInvers
*
* INPUT
*
*   m - matrix to invert
*   r - inverted matrix
*   
* OUTPUT
*
*   r
*   
* RETURNS
*   
* AUTHOR
*
*   Dieter Bayer
*   
* DESCRIPTION
*
*   Invert a 4x4 matrix.
*
* CHANGES
*
*   May 1994 : Creation.
*
******************************************************************************/

void MInvers(MATRIX r, const MATRIX  m)
{
	DBL d00, d01, d02, d03;
	DBL d10, d11, d12, d13;
	DBL d20, d21, d22, d23;
	DBL d30, d31, d32, d33;
	DBL m00, m01, m02, m03;
	DBL m10, m11, m12, m13;
	DBL m20, m21, m22, m23;
	DBL m30, m31, m32, m33;
	DBL D;

	m00 = m[0][0];  m01 = m[0][1];  m02 = m[0][2];  m03 = m[0][3];
	m10 = m[1][0];  m11 = m[1][1];  m12 = m[1][2];  m13 = m[1][3];
	m20 = m[2][0];  m21 = m[2][1];  m22 = m[2][2];  m23 = m[2][3];
	m30 = m[3][0];  m31 = m[3][1];  m32 = m[3][2];  m33 = m[3][3];

	d00 = m11*m22*m33 + m12*m23*m31 + m13*m21*m32 - m31*m22*m13 - m32*m23*m11 - m33*m21*m12;
	d01 = m10*m22*m33 + m12*m23*m30 + m13*m20*m32 - m30*m22*m13 - m32*m23*m10 - m33*m20*m12;
	d02 = m10*m21*m33 + m11*m23*m30 + m13*m20*m31 - m30*m21*m13 - m31*m23*m10 - m33*m20*m11;
	d03 = m10*m21*m32 + m11*m22*m30 + m12*m20*m31 - m30*m21*m12 - m31*m22*m10 - m32*m20*m11;

	d10 = m01*m22*m33 + m02*m23*m31 + m03*m21*m32 - m31*m22*m03 - m32*m23*m01 - m33*m21*m02;
	d11 = m00*m22*m33 + m02*m23*m30 + m03*m20*m32 - m30*m22*m03 - m32*m23*m00 - m33*m20*m02;
	d12 = m00*m21*m33 + m01*m23*m30 + m03*m20*m31 - m30*m21*m03 - m31*m23*m00 - m33*m20*m01;
	d13 = m00*m21*m32 + m01*m22*m30 + m02*m20*m31 - m30*m21*m02 - m31*m22*m00 - m32*m20*m01;

	d20 = m01*m12*m33 + m02*m13*m31 + m03*m11*m32 - m31*m12*m03 - m32*m13*m01 - m33*m11*m02;
	d21 = m00*m12*m33 + m02*m13*m30 + m03*m10*m32 - m30*m12*m03 - m32*m13*m00 - m33*m10*m02;
	d22 = m00*m11*m33 + m01*m13*m30 + m03*m10*m31 - m30*m11*m03 - m31*m13*m00 - m33*m10*m01;
	d23 = m00*m11*m32 + m01*m12*m30 + m02*m10*m31 - m30*m11*m02 - m31*m12*m00 - m32*m10*m01;

	d30 = m01*m12*m23 + m02*m13*m21 + m03*m11*m22 - m21*m12*m03 - m22*m13*m01 - m23*m11*m02;
	d31 = m00*m12*m23 + m02*m13*m20 + m03*m10*m22 - m20*m12*m03 - m22*m13*m00 - m23*m10*m02;
	d32 = m00*m11*m23 + m01*m13*m20 + m03*m10*m21 - m20*m11*m03 - m21*m13*m00 - m23*m10*m01;
	d33 = m00*m11*m22 + m01*m12*m20 + m02*m10*m21 - m20*m11*m02 - m21*m12*m00 - m22*m10*m01;

	D = m00*d00 - m01*d01 + m02*d02 - m03*d03;

	if (D == 0.0)
	{
		throw POV_EXCEPTION_STRING("Singular matrix in MInvers.");
	}

	r[0][0] =  d00/D; r[0][1] = -d10/D;  r[0][2] =  d20/D; r[0][3] = -d30/D;
	r[1][0] = -d01/D; r[1][1] =  d11/D;  r[1][2] = -d21/D; r[1][3] =  d31/D;
	r[2][0] =  d02/D; r[2][1] = -d12/D;  r[2][2] =  d22/D; r[2][3] = -d32/D;
	r[3][0] = -d03/D; r[3][1] =  d13/D;  r[3][2] = -d23/D; r[3][3] =  d33/D;
}

UV_VECT *Create_UV_Vect ()
{
	UV_VECT *New;

	New = reinterpret_cast<UV_VECT *>(POV_MALLOC(sizeof (UV_VECT), "uv vector"));

	(*New)[0]= 0.0;
	(*New)[1]= 0.0;

	return (New);
}

VECTOR_4D *Create_Vector_4D ()
{
	VECTOR_4D *New;

	New = reinterpret_cast<VECTOR_4D *>(POV_MALLOC(sizeof (VECTOR_4D), "4d vector"));

	(*New)[0]= 0.0;
	(*New)[1]= 0.0;
	(*New)[2]= 0.0;
	(*New)[3]= 0.0;

	return (New);
}

/*****************************************************************************
*
* FUNCTION
*
*   MTransUVPoint
*
* INPUT
*
*   p - point to tranform
*   m - 3x3 matrix
*   
* OUTPUT
*
*   t
*   
* RETURNS
*   
* AUTHOR
*
*   Michael Hough
*   
* DESCRIPTION
*
* CHANGES
*
*   June 2000 : Creation.
*
******************************************************************************/

void MTransUVPoint(const DBL p[2], const DBL m[3][3], DBL t[2])
{
	DBL w;

	t[0] = p[0]*m[0][0] + p[1]*m[1][0] + m[2][0];
	t[1] = p[0]*m[0][1] + p[1]*m[1][1] + m[2][1];

	w = p[0]*m[0][2] + p[1]*m[1][2] + m[2][2];

	t[0] /= w;
	t[1] /= w;
}


/*****************************************************************************
*
* FUNCTION
*
*   MSquare_To_Quad
*
* INPUT
*
*   st - four texture coordinates
*   
* OUTPUT
*
*   sq - 3x3 matrix
*   
* RETURNS
*   
* AUTHOR
*
*   Michael Hough
* 
* DESCRIPTION
*
*   Find the forward mapping matrix - unit square -> general quadrilateral
*
*   formula from: "Fundamentals of Texture Mapping and Image Warping", Paul Heckbert
*
* CHANGES
*
*   June 2000 : Creation.
*
******************************************************************************/

void MSquareQuad(const UV_VECT st[4], DBL sq[3][3])
{
	DBL sx, sy, dx1, dx2, dy1, dy2, det;

	dx1 = st[1][0] - st[2][0];
	dx2 = st[3][0] - st[2][0];
	dy1 = st[1][1] - st[2][1];
	dy2 = st[3][1] - st[2][1];

	sx = st[0][0] - st[1][0] + st[2][0] - st[3][0];
	sy = st[0][1] - st[1][1] + st[2][1] - st[3][1];

	det = Det(dx1, dx2, dy1, dy2);

	if(det == 0.0)
	{
		throw POV_EXCEPTION_STRING("Equal uv-vectors detected, division by zero!");
	}

	sq[0][2] = Det(sx, dx2, sy, dy2)/det;
	sq[1][2] = Det(dx1, sx, dy1, sy)/det;
	sq[2][2] = 1;
	sq[0][0] = st[1][0] - st[0][0] + sq[0][2]*st[1][0];
	sq[1][0] = st[3][0] - st[0][0] + sq[1][2]*st[3][0];
	sq[2][0] = st[0][0];
	sq[0][1] = st[1][1] - st[0][1] + sq[0][2]*st[1][1];
	sq[1][1] = st[3][1] - st[0][1] + sq[1][2]*st[3][1];
	sq[2][1] = st[0][1];
}

}
