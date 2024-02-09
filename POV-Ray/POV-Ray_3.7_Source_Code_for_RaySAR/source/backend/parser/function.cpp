/*******************************************************************************
 * function.cpp
 *
 * This module implements the the function type used by iso surfaces and
 * the function pattern.
 *
 * This module is based on code by D. Skarda, T. Bily and R. Suzuki.
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
 * $File: //depot/public/povray/3.x/source/backend/parser/function.cpp $
 * $Revision: #1 $
 * $Change: 6069 $
 * $DateTime: 2013/11/06 11:59:40 $
 * $Author: chrisc $
 *******************************************************************************/

#include <limits.h>

// frame.h must always be the first POV file included (pulls in platform config)
#include "backend/frame.h"
#include "backend/parser/parse.h"
#include "backend/vm/fnpovfpu.h"
#include "backend/math/splines.h"
#include "backend/texture/pigment.h"

// this must be the last file included
#include "base/povdebug.h"

namespace pov
{

/*****************************************************************************
 *
 * FUNCTION
 *
 *   Destroy_Function
 *
 * INPUT
 *
 * OUTPUT
 *
 * RETURNS
 *
 * AUTHOR
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

void Parser::Destroy_Function(FUNCTION_PTR Function)
{
	if(Function != NULL)
	{
		sceneData->functionVM->RemoveFunction(*Function);
		POV_FREE(Function);
	}
}

void Parser::Destroy_Function(FunctionVM *functionVM, FUNCTION_PTR Function)
{
	if(Function != NULL)
	{
		functionVM->RemoveFunction(*Function);
		POV_FREE(Function);
	}
}


/*****************************************************************************
 *
 * FUNCTION
 *
 *   Copy_Function
 *
 * INPUT
 *
 * OUTPUT
 *
 * RETURNS
 *
 * AUTHOR
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

FUNCTION_PTR Parser::Copy_Function(FUNCTION_PTR Function)
{
	FUNCTION_PTR ptr = (FUNCTION_PTR)POV_MALLOC(sizeof(FUNCTION), "Function ID");

	(void)sceneData->functionVM->GetFunctionAndReference(*Function); // increase the reference count
	*ptr = *Function;

	return ptr;
}

FUNCTION_PTR Parser::Copy_Function(FunctionVM *functionVM, FUNCTION_PTR Function)
{
	FUNCTION_PTR ptr = (FUNCTION_PTR)POV_MALLOC(sizeof(FUNCTION), "Function ID");

	functionVM->GetFunctionAndReference(*Function); // increase the reference count
	*ptr = *Function;

	return ptr;
}


/*****************************************************************************
*
* FUNCTION
*
*   Parse_Function
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*
*   FUNCTION - parsed and compiled function reference number
*
* AUTHOR
*
*   Thorsten Froehlich
*   
* DESCRIPTION
*
*   Parse and compile a function and add it to the global function table.
*
* CHANGES
*
*   -
*
******************************************************************************/

FUNCTION_PTR Parser::Parse_Function(void)
{
	FUNCTION_PTR ptr = (FUNCTION_PTR)POV_MALLOC(sizeof(FUNCTION), "Function ID");
	ExprNode *expression = NULL;
	FunctionCode function;

	Parse_Begin();

	FNCode f(this, &function, false, NULL);

	expression = FNSyntax_ParseExpression();
	f.Compile(expression);
	FNSyntax_DeleteExpression(expression);

	Parse_End();

	*ptr = sceneData->functionVM->AddFunction(&function);

	return ptr;
}


/*****************************************************************************
*
* FUNCTION
*
*   Parse_FunctionContent
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*
*   FUNCTION - parsed and compiled function reference number
*
* AUTHOR
*
*   Thorsten Froehlich
*   
* DESCRIPTION
*
*   Parse and compile a function and add it to the global function table.
*
* CHANGES
*
*   -
*
******************************************************************************/

FUNCTION_PTR Parser::Parse_FunctionContent(void)
{
	FUNCTION_PTR ptr = (FUNCTION_PTR)POV_MALLOC(sizeof(FUNCTION), "Function ID");
	ExprNode *expression = NULL;
	FunctionCode function;

	FNCode f(this, &function, false, NULL);

	expression = FNSyntax_ParseExpression();
	f.Compile(expression);
	FNSyntax_DeleteExpression(expression);

	*ptr = sceneData->functionVM->AddFunction(&function);

	return ptr;
}


/*****************************************************************************
*
* FUNCTION
*
*   Parse_DeclareFunction
*
* INPUT
*   
* OUTPUT
*   
* RETURNS
*
*   FUNCTION - parsed and compiled function reference number
*
* AUTHOR
*
*   Thorsten Froehlich
*   
* DESCRIPTION
*
*   Parse and compile a function and add it to the global function table.
*   Additionally, this function takes an optional list of parameter names.
*
* CHANGES
*
*   -
*
******************************************************************************/

FUNCTION_PTR Parser::Parse_DeclareFunction(int *token_id, const char *fn_name, bool is_local)
{
	FUNCTION_PTR ptr = (FUNCTION_PTR)POV_MALLOC(sizeof(FUNCTION), "Function ID");
	ExprNode *expression = NULL;
	FunctionCode function;

	// default type is float function
	*token_id = FUNCT_ID_TOKEN;

	FNCode f(this, &function, is_local, fn_name);
	f.Parameter();

	Parse_Begin();

	Get_Token();
	if(Token.Token_Id == INTERNAL_TOKEN)
	{
		GET(LEFT_PAREN_TOKEN);

		Get_Token();
		if(Token.Function_Id != FLOAT_TOKEN)
			Expectation_Error("internal function identifier");
		expression = FNSyntax_GetTrapExpression((unsigned int)(Token.Token_Float));

		function.flags = FN_INLINE_FLAG;

		GET(RIGHT_PAREN_TOKEN);
	}
	else if(Token.Token_Id == TRANSFORM_TOKEN)
	{
		if(function.parameter_cnt != 0)
			Error("Function parameters for transform functions are not allowed.");

		expression = FNSyntax_GetTrapExpression(1); // 1 refers to POVFPU_TrapSTable[1] = f_transform [trf]

		function.private_copy_method = (FNCODE_PRIVATE_COPY_METHOD)Copy_Transform;
		function.private_destroy_method = (FNCODE_PRIVATE_DESTROY_METHOD)Destroy_Transform;

		function.private_data = (void *)Parse_Transform_Block();

		function.return_size = 3; // returns a 3d vector!!!

		// function type is vector function
		*token_id = VECTFUNCT_ID_TOKEN;
	}
	else if(Token.Token_Id == SPLINE_TOKEN)
	{
		if(function.parameter_cnt != 0)
			Error("Function parameters for spline functions are not allowed.");

		Experimental_Flag |= EF_SPLINE;

		expression = FNSyntax_GetTrapExpression(2); // 2 refers to POVFPU_TrapSTable[2] = f_spline [trf]

		function.private_copy_method = (FNCODE_PRIVATE_COPY_METHOD)Copy_Spline;
		function.private_destroy_method = (FNCODE_PRIVATE_DESTROY_METHOD)Destroy_Spline;

		Parse_Begin();
		function.private_data = (void *)Parse_Spline();
		Parse_End();

		function.return_size = ((SPLINE *)(function.private_data))->Terms; // returns a 2d, 3d, 4d or 5d vector!!!

		// function type is vector function
		*token_id = VECTFUNCT_ID_TOKEN;
	}
	else if(Token.Token_Id == PIGMENT_TOKEN)
	{
		if(function.parameter_cnt != 0)
			Error("Function parameters for pigment functions are not allowed.");

		expression = FNSyntax_GetTrapExpression(0); // 0 refers to POVFPU_TrapSTable[0] = f_pigment [trf]

		function.private_copy_method = (FNCODE_PRIVATE_COPY_METHOD)Copy_Pigment;
		function.private_destroy_method = (FNCODE_PRIVATE_DESTROY_METHOD)Destroy_Pigment;

		Parse_Begin();
		function.private_data = (void *)Create_Pigment();
		Parse_Pigment((PIGMENT **)(&function.private_data));
		Parse_End();
		Post_Pigment((PIGMENT *)(function.private_data));

		function.return_size = 5; // returns a color!!!

		// function type is vector function
		*token_id = VECTFUNCT_ID_TOKEN;
	}
	else if(Token.Token_Id == PATTERN_TOKEN)
	{
		if(function.parameter_cnt != 0)
			Error("Function parameters for pattern functions are not allowed.");

		expression = FNSyntax_GetTrapExpression(77); // 77 refers to POVFPU_TrapTable[77] = f_pattern [trf]

		function.private_copy_method = (FNCODE_PRIVATE_COPY_METHOD)Copy_Pigment;
		function.private_destroy_method = (FNCODE_PRIVATE_DESTROY_METHOD)Destroy_Pigment;

		Parse_Begin();
		function.private_data = (void *)Create_Pigment(); // Yes, this is a pigment! [trf]
		Parse_PatternFunction((TPATTERN *)(function.private_data));
		Parse_End();
		Post_Pigment((PIGMENT *)(function.private_data));
	}
	else if(Token.Token_Id == STRING_LITERAL_TOKEN)
	{
		f.SetFlag(2, Token.Token_String);
		Get_Token();
		if(Token.Token_Id == COMMA_TOKEN)
		{
			Get_Token();
			if(Token.Token_Id != STRING_LITERAL_TOKEN)
				Expectation_Error("valid function expression");
			f.SetFlag(1, Token.Token_String);
		}
		else
		{
			Unget_Token();
			expression = FNSyntax_ParseExpression();
		}
	}
	else
	{
		Unget_Token();
		expression = FNSyntax_ParseExpression();
	}

	f.Compile(expression);
	FNSyntax_DeleteExpression(expression);

	Parse_End();

	*ptr = sceneData->functionVM->AddFunction(&function);

	return ptr;
}

}
