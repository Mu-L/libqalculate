/*
    Qalculate (library)

    Copyright (C) 2003-2007, 2008, 2016, 2018  Hanna Knutsson (hanna.knutsson@protonmail.com)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
*/

#include "support.h"

#include "BuiltinFunctions.h"
#include "util.h"
#include "MathStructure.h"
#include "Number.h"
#include "Calculator.h"
#include "Variable.h"
#include "Unit.h"

#include <sstream>
#include <time.h>
#include <limits>
#include <algorithm>

#include "MathStructure-support.h"

using std::string;
using std::cout;
using std::vector;
using std::endl;

VectorFunction::VectorFunction() : MathFunction("vector", -1) {
}
int VectorFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct = vargs;
	mstruct.setType(STRUCT_VECTOR);
	return 1;
}
MatrixFunction::MatrixFunction() : MathFunction("matrix", 3) {
	Argument *arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	arg->setHandleVector(false);
	setArgumentDefinition(1, arg);
	arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	arg->setHandleVector(false);
	setArgumentDefinition(2, arg);
	setArgumentDefinition(3, new VectorArgument());
}
int MatrixFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	size_t rows = (size_t) vargs[0].number().uintValue();
	size_t columns = (size_t) vargs[1].number().uintValue();
	if((rows > 1000 || columns > 1000) && vargs[0].number() * vargs[1].number() > Number(1, 1, 6)) return 0;
	mstruct.clearMatrix(); mstruct.resizeMatrix(rows, columns, m_zero);
	if(mstruct.rows() < rows || mstruct.columns() < columns) return 0;
	size_t r = 1, c = 1;
	for(size_t i = 0; i < vargs[2].size(); i++) {
		if(r > rows || c > columns) {
			CALCULATOR->error(false, _("Too many elements (%s) for the dimensions (%sx%s) of the matrix."), i2s(vargs[2].size()).c_str(), i2s(rows).c_str(), i2s(columns).c_str(), NULL);
			break;
		}
		mstruct[r - 1][c - 1] = vargs[2][i];
		if(c == columns) {
			c = 1;
			r++;
		} else {
			c++;
		}
	}
	return 1;
}
RankFunction::RankFunction() : MathFunction("rank", 1, 2) {
	setArgumentDefinition(1, new VectorArgument(""));
	setArgumentDefinition(2, new BooleanArgument(""));
	setDefaultValue(2, "1");
}
int RankFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	if(vargs[0].isMatrix()) {
		MathStructure mvector;
		mvector.clearVector();
		size_t rows = vargs[0].size(), cols = vargs[0][0].size();
		for(size_t i = 0; i < rows; i++) {
			for(size_t i2 = 0; i2 < cols; i2++) {
				if(CALCULATOR->aborted()) return 0;
				mvector.addChild(vargs[0][i][i2]);
			}
		}
		if(!mvector.rankVector(vargs[1].number().getBoolean())) return 0;
		mstruct.clearMatrix();
		mstruct.resizeMatrix(rows, cols, m_zero);
		if(mstruct.rows() < rows || mstruct.columns() < cols) return 0;
		for(size_t i = 0; i < mvector.size(); i++) {
			if(CALCULATOR->aborted()) return 0;
			mstruct[i / cols][i % cols] = mvector[i];
		}
		return 1;
	}
	mstruct = vargs[0];
	return mstruct.rankVector(vargs[1].number().getBoolean());
}
SortFunction::SortFunction() : MathFunction("sort", 1, 2) {
	setArgumentDefinition(1, new VectorArgument(""));
	setArgumentDefinition(2, new BooleanArgument(""));
	setDefaultValue(2, "1");
}
int SortFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct = vargs[0];
	return mstruct.sortVector(vargs[1].number().getBoolean());
}
MergeVectorsFunction::MergeVectorsFunction() : MathFunction("mergevectors", 1, -1) {
	setArgumentDefinition(1, new VectorArgument(""));
	setArgumentDefinition(2, new VectorArgument(""));
}
int MergeVectorsFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct.clearVector();
	for(size_t i = 0; i < vargs.size(); i++) {
		if(vargs[i].isVector()) {
			for(size_t i2 = 0; i2 < vargs[i].size(); i2++) {
				if(CALCULATOR->aborted()) return 0;
				mstruct.addChild(vargs[i][i2]);
			}
		} else {
			if(CALCULATOR->aborted()) return 0;
			mstruct.addChild(vargs[i]);
		}
	}
	return 1;
}
FlipFunction::FlipFunction() : MathFunction("flip", 1, 2) {
	setArgumentDefinition(1, new MatrixArgument());
	IntegerArgument *iarg = new IntegerArgument();
	iarg->setMin(&nr_zero);
	iarg->setMax(&nr_two);
	setArgumentDefinition(2, iarg);
	setDefaultValue(2, "0");
}
int FlipFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct = vargs[0];
	int dimension = vargs[1].number().intValue();
	if(dimension == 0 || dimension == 1) {
		mstruct.flipVector();
	}
	if(dimension == 0 || dimension == 2) {
		for(size_t i = 0; i < mstruct.size(); i++) {
			mstruct[i].flipVector();
		}
	}
	return 1;
}
VertCatFunction::VertCatFunction() : MathFunction("vertcat", 1, -1) {
	setArgumentDefinition(1, new MatrixArgument(""));
	setArgumentDefinition(2, new MatrixArgument(""));
}
int VertCatFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct = vargs[0];
	for(size_t i = 1; i < vargs.size(); i++) {
		if(vargs[i].columns() != mstruct.columns()) {
			CALCULATOR->error(true, _("Vertical concatenation requires equal number of columns."), NULL);
			if(i > 1) {
				mstruct.transform(this);
				for(; i < vargs.size(); i++) mstruct.addChild(vargs[i]);
				return 1;
			}
			return 0;
		}
		for(size_t i2 = 0; i2 < vargs[i].size(); i2++) {
			if(CALCULATOR->aborted()) return 0;
			mstruct.addChild(vargs[i][i2]);
		}
	}
	return 1;
}
HorzCatFunction::HorzCatFunction() : MathFunction("horzcat", 1, -1) {
	setArgumentDefinition(1, new MatrixArgument(""));
	setArgumentDefinition(2, new MatrixArgument(""));
}
bool HorzCatFunction::representsScalar(const MathStructure &vargs) const {return false;}
bool HorzCatFunction::representsNonMatrix(const MathStructure &vargs) const {
	for(size_t i = 0; i < vargs.size(); i++) {
		if(!vargs[i].representsNonMatrix()) return false;
	}
	return true;
}
int HorzCatFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct = vargs[0];
	for(size_t i = 1; i < vargs.size(); i++) {
		if(vargs[i].rows() != mstruct.rows()) {
			CALCULATOR->error(true, _("Horizontal concatenation requires equal number of rows."), NULL);
			if(i > 1) {
				mstruct.transform(this);
				for(; i < vargs.size(); i++) mstruct.addChild(vargs[i]);
				return 1;
			}
			return 0;
		}
		for(size_t i2 = 0; i2 < vargs[i].size(); i2++) {
			for(size_t i3 = 0; i3 < vargs[i][i2].size(); i3++) {
				if(CALCULATOR->aborted()) return 0;
				mstruct[i2].addChild(vargs[i][i2][i3]);
			}
		}
	}
	return 1;
}
MatrixToVectorFunction::MatrixToVectorFunction() : MathFunction("matrix2vector", 1) {
	setArgumentDefinition(1, new MatrixArgument());
}
int MatrixToVectorFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	vargs[0].matrixToVector(mstruct);
	return 1;
}
RowFunction::RowFunction() : MathFunction("row", 2) {
	setArgumentDefinition(1, new MatrixArgument());
	setArgumentDefinition(2, new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE));
}
int RowFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	size_t row = (size_t) vargs[1].number().uintValue();
	if(row > vargs[0].rows()) {
		CALCULATOR->error(true, _("Row %s does not exist in matrix."), format_and_print(vargs[1]).c_str(), NULL);
		return 0;
	}
	vargs[0].rowToVector(row, mstruct);
	return 1;
}
ColumnFunction::ColumnFunction() : MathFunction("column", 2) {
	setArgumentDefinition(1, new MatrixArgument());
	setArgumentDefinition(2, new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE));
}
int ColumnFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	size_t col = (size_t) vargs[1].number().uintValue();
	if(col > vargs[0].columns()) {
		CALCULATOR->error(true, _("Column %s does not exist in matrix."), format_and_print(vargs[1]).c_str(), NULL);
		return 0;
	}
	vargs[0].columnToVector(col, mstruct);
	return 1;
}
RowsFunction::RowsFunction() : MathFunction("rows", 1) {
	setArgumentDefinition(1, new MatrixArgument(""));
}
int RowsFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct.set((long int) vargs[0].rows(), 1L, 0L);
	return 1;
}
ColumnsFunction::ColumnsFunction() : MathFunction("columns", 1) {
	setArgumentDefinition(1, new MatrixArgument(""));
}
int ColumnsFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct.set((long int) vargs[0].columns(), 1L, 0L);
	return 1;
}
ElementsFunction::ElementsFunction() : MathFunction("elements", 1) {
	setArgumentDefinition(1, new MatrixArgument("", false));
}
int ElementsFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct = vargs[0];
	mstruct.eval(eo);
	if(!mstruct.isMatrix()) {
		if(mstruct.isVector() && (mstruct.size() == 0 || mstruct[0].representsScalar())) {
			mstruct.set((long int) mstruct.size(), 1L, 0L);
			return 1;
		} else if(mstruct.representsScalar()) {
			mstruct = m_one;
			return 1;
		} else if(eo.approximation == APPROXIMATION_EXACT || eo.approximation == APPROXIMATION_EXACT_VARIABLES) {
			EvaluationOptions eo2 = eo;
			eo2.approximation = APPROXIMATION_APPROXIMATE;
			MathStructure m2(vargs[0]);
			CALCULATOR->beginTemporaryStopMessages();
			m2.eval(eo2);
			if(CALCULATOR->endTemporaryStopMessages()) return -1;
			if(m2.isMatrix()) {
				mstruct.set((long int) (m2.rows() * m2.columns()), 1L, 0L);
				return 1;
			} else if(m2.isVector() && (m2.size() == 0 || m2[0].representsScalar())) {
				mstruct.set((long int) m2.size(), 1L, 0L);
				return 1;
			} else if(m2.representsScalar()) {
				mstruct = m_one;
				return 1;
			}
		}
		return -1;
	}
	mstruct.set((long int) (mstruct.rows() * mstruct.columns()), 1L, 0L);
	return 1;
}
ElementFunction::ElementFunction() : MathFunction("element", 2, 3) {
	setArgumentDefinition(1, new MatrixArgument(""));
	setArgumentDefinition(2, new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE));
	setArgumentDefinition(3, new IntegerArgument("", ARGUMENT_MIN_MAX_NONE, true, true, INTEGER_TYPE_SIZE));
	setDefaultValue(3, "0");
}
bool ElementFunction::representsScalar(const MathStructure &vargs) const {
	if(vargs.size() >= 2 && vargs[0].isMatrix() && vargs[1].isInteger() && vargs[1].number().isPositive()) {
		if(vargs.size() == 2 || vargs[2].isZero()) {
			if(vargs[0].size() == 1 && vargs[1].number() <= vargs[0][0].size()) return vargs[0][0][vargs[1].number().uintValue() - 1].representsScalar();
			if(vargs[1].number() <= vargs[0].size() && vargs[0][0].size() == 1) return vargs[0][vargs[1].number().uintValue() - 1][0].representsScalar();
		} else if(vargs[1].number() <= vargs[0].size() && vargs[2].isInteger() && vargs[2].number().isPositive() && vargs[2].number() <= vargs[0][0].size()) {
			return vargs[0][vargs[1].number().uintValue() - 1][vargs[2].number().uintValue() - 1].representsScalar();
		}
	}
	return false;
}
bool ElementFunction::representsNonMatrix(const MathStructure &vargs) const {
	if(vargs.size() >= 2 && vargs[0].isMatrix() && vargs[1].isInteger() && vargs[1].number().isPositive()) {
		if(vargs.size() == 2 || vargs[2].isZero()) {
			if(vargs[0].size() == 1 && vargs[1].number() <= vargs[0][0].size()) return vargs[0][0][vargs[1].number().uintValue() - 1].representsNonMatrix();
			if(vargs[1].number() <= vargs[0].size()) {
				if(vargs[0][0].size() == 1) return vargs[0][vargs[1].number().uintValue() - 1][0].representsNonMatrix();
				else return vargs[0][vargs[1].number().uintValue() - 1].representsNonMatrix();
			}
		} else if(vargs[1].number() <= vargs[0].size() && vargs[2].isInteger() && vargs[2].number().isPositive() && vargs[2].number() <= vargs[0][0].size()) {
			return vargs[0][vargs[1].number().uintValue() - 1][vargs[2].number().uintValue() - 1].representsNonMatrix();
		}
	}
	return false;
}
int ElementFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	size_t row = (size_t) vargs[1].number().uintValue();
	size_t col = (size_t) vargs[2].number().uintValue();
	if(col == 0) {
		if(vargs[0].size() == 1 && vargs[0][0].size() >= row) {
			mstruct = vargs[0][0][row - 1];
			return 1;
		}
		if(vargs[0].size() >= row) {
			if(vargs[0][0].size() == 1) mstruct = vargs[0][row - 1][0];
			else mstruct = vargs[0][row - 1];
			return 1;
		}
		CALCULATOR->error(true, _("Element %s does not exist in vector."), format_and_print(vargs[1]).c_str(), NULL);
		return 0;
	}
	if(col > vargs[0].columns()) {
		CALCULATOR->error(true, _("Column %s does not exist in matrix."), format_and_print(vargs[2]).c_str(), NULL);
		return 0;
	}
	if(row > vargs[0].rows()) {
		CALCULATOR->error(true, _("Row %s does not exist in matrix."), format_and_print(vargs[1]).c_str(), NULL);
		return 0;
	}
	const MathStructure *em = vargs[0].getElement(row, col);
	if(em) mstruct = *em;
	else return 0;
	return 1;
}
DimensionFunction::DimensionFunction() : MathFunction("dimension", 1) {
	setArgumentDefinition(1, new VectorArgument(""));
}
int DimensionFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct.set((long int) vargs[0].countChildren(), 1L, 0L);
	return 1;
}
ComponentFunction::ComponentFunction() : MathFunction("component", 2) {
	setArgumentDefinition(1, new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE));
	setArgumentDefinition(2, new VectorArgument(""));
}
bool ComponentFunction::representsScalar(const MathStructure &vargs) const {
	if(vargs.size() >= 2 && vargs[0].isVector() && vargs[1].isInteger() && vargs[1].number().isPositive() && vargs[1].number() <= vargs[0].size()) {
		return vargs[0][vargs[1].number().uintValue() - 1].representsScalar();
	}
	return false;
}
bool ComponentFunction::representsNonMatrix(const MathStructure &vargs) const {
	if(vargs.size() >= 2 && vargs[0].isVector() && vargs[1].isInteger() && vargs[1].number().isPositive() && vargs[1].number() <= vargs[0].size()) {
		return vargs[0][vargs[1].number().uintValue() - 1].representsNonMatrix();
	}
	return false;
}
int ComponentFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	size_t i = (size_t) vargs[0].number().uintValue();
	if(i > vargs[1].countChildren()) {
		CALCULATOR->error(true, _("Element %s does not exist in vector."), format_and_print(vargs[0]).c_str(), NULL);
		return 0;
	}
	mstruct = *vargs[1].getChild(i);
	return 1;
}
LimitsFunction::LimitsFunction() : MathFunction("limits", 3) {
	setArgumentDefinition(1, new VectorArgument(""));
	Argument *arg = new IntegerArgument("", ARGUMENT_MIN_MAX_NONE, true, true, INTEGER_TYPE_SINT);
	arg->setHandleVector(false);
	setArgumentDefinition(2, arg);
	arg = new IntegerArgument("", ARGUMENT_MIN_MAX_NONE, true, true, INTEGER_TYPE_SINT);
	arg->setHandleVector(false);
	setArgumentDefinition(3, arg);
}
int LimitsFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	vargs[0].getRange(vargs[1].number().intValue(), vargs[2].number().intValue(), mstruct);
	return 1;
}
AreaFunction::AreaFunction() : MathFunction("area", 5) {
	setArgumentDefinition(1, new MatrixArgument(""));
	Argument *arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	arg->setHandleVector(false);
	setArgumentDefinition(2, arg);
	arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	arg->setHandleVector(false);
	setArgumentDefinition(3, arg);
	arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	arg->setHandleVector(false);
	setArgumentDefinition(4, arg);
	arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	arg->setHandleVector(false);
	setArgumentDefinition(5, arg);
}
int AreaFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	vargs[0].getArea(vargs[1].number().uintValue(), vargs[2].number().uintValue(), vargs[3].number().uintValue(), vargs[4].number().uintValue(), mstruct);
	if(vargs[0].isUndefined()) return 0;
	return 1;
}
TransposeFunction::TransposeFunction() : MathFunction("transpose", 1) {
	setArgumentDefinition(1, new MatrixArgument());
}
int TransposeFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	mstruct = vargs[0];
	return mstruct.transposeMatrix();
}
IdentityMatrixFunction::IdentityMatrixFunction() : MathFunction("identity", 1) {
	ArgumentSet *arg = new ArgumentSet();
	IntegerArgument *iarg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE);
	Number nr(1, 1, 7);
	iarg->setMax(&nr);
	arg->addArgument(iarg);
	MatrixArgument *marg = new MatrixArgument();
	marg->setSquareDemanded(true);
	arg->addArgument(marg);
	setArgumentDefinition(1, arg);
}
int IdentityMatrixFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	if(vargs[0].isMatrix()) {
		if(vargs[0].rows() != vargs[0].columns()) {
			return 0;
		}
		mstruct.setToIdentityMatrix(vargs[0].size());
	} else {
		mstruct.setToIdentityMatrix((size_t) vargs[0].number().uintValue());
	}
	if(mstruct.isUndefined()) return 0;
	return 1;
}
DeterminantFunction::DeterminantFunction() : MathFunction("det", 1) {
	MatrixArgument *marg = new MatrixArgument();
	marg->setSquareDemanded(true);
	setArgumentDefinition(1, marg);
}
int DeterminantFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	vargs[0].determinant(mstruct, eo);
	return !mstruct.isUndefined();
}
PermanentFunction::PermanentFunction() : MathFunction("permanent", 1) {
	MatrixArgument *marg = new MatrixArgument();
	marg->setSquareDemanded(true);
	setArgumentDefinition(1, marg);
}
int PermanentFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	vargs[0].permanent(mstruct, eo);
	return !mstruct.isUndefined();
}
CofactorFunction::CofactorFunction() : MathFunction("cofactor", 3) {
	MatrixArgument *marg = new MatrixArgument();
	marg->setSquareDemanded(true);
	setArgumentDefinition(1, marg);
	setArgumentDefinition(2, new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE));
	setArgumentDefinition(3, new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SIZE));
}
int CofactorFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	vargs[0].cofactor((size_t) vargs[1].number().uintValue(), (size_t) vargs[2].number().uintValue(), mstruct, eo);
	return !mstruct.isUndefined();
}
AdjointFunction::AdjointFunction() : MathFunction("adj", 1) {
	MatrixArgument *marg = new MatrixArgument();
	marg->setSquareDemanded(true);
	setArgumentDefinition(1, marg);
}
int AdjointFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct = vargs[0];
	if(!mstruct.adjointMatrix(eo)) return 0;
	return !mstruct.isUndefined();
}
InverseFunction::InverseFunction() : MathFunction("inv", 1) {
	MatrixArgument *marg = new MatrixArgument();
	marg->setTests(false);
	marg->setSquareDemanded(true);
	setArgumentDefinition(1, marg);
}
int InverseFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct = vargs[0];
	if(mstruct.representsScalar()) {mstruct.inverse(); return 1;}
	if(!mstruct.isMatrix()) {
		mstruct.eval(eo);
		if(mstruct.representsScalar()) {mstruct.inverse(); return 1;}
		else if(!mstruct.isVector()) return -1;
	}
	if(!mstruct.isMatrix() || !mstruct.matrixIsSquare()) {
		Argument *arg = getArgumentDefinition(1);
		arg->setTests(true);
		arg->test(mstruct, 1, this, eo);
		arg->setTests(false);
		return -1;
	}
	return mstruct.invertMatrix(eo);
}
MagnitudeFunction::MagnitudeFunction() : MathFunction("magnitude", 1) {
	setArgumentDefinition(1, new VectorArgument(""));
}
int MagnitudeFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	if(vargs[0].isMatrix()) return 0;
	mstruct = vargs[0];
	if(!mstruct.representsNonMatrix()) {
		CALCULATOR->beginTemporaryStopMessages();
		mstruct.eval(eo);
		CALCULATOR->endTemporaryStopMessages();
		if(mstruct.isMatrix()) return 0;
		mstruct = vargs[0];
	}
	if(mstruct.size() == 1) {
		mstruct.setType(STRUCT_FUNCTION);
		mstruct.setFunctionId(FUNCTION_ID_ABS);
	} else {
		for(size_t i = 0; i < mstruct.size(); i++) {
			if(!mstruct[i].representsReal(true)) mstruct[i].transformById(FUNCTION_ID_ABS);
			mstruct[i] ^= 2;
		}
		if(mstruct.size() == 0) mstruct.clear(true);
		else if(mstruct.size() == 1) mstruct.setToChild(1, true);
		else mstruct.setType(STRUCT_ADDITION);
		mstruct ^= nr_half;
	}
	return 1;
}
NormFunction::NormFunction() : MathFunction("norm", 1, 2) {
	setArgumentDefinition(1, new VectorArgument(""));
	setArgumentDefinition(2, new IntegerArgument("", ARGUMENT_MIN_MAX_NONE, false, false));
	setDefaultValue(2, "2");
}
int NormFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	if(vargs[0].isMatrix()) return 0;
	mstruct = vargs[0];
	if(!mstruct.representsNonMatrix()) {
		CALCULATOR->beginTemporaryStopMessages();
		mstruct.eval(eo);
		CALCULATOR->endTemporaryStopMessages();
		if(mstruct.isMatrix()) return 0;
		mstruct = vargs[0];
	}
	if(mstruct.size() == 1) {
		mstruct.setType(STRUCT_FUNCTION);
		mstruct.setFunctionId(FUNCTION_ID_ABS);
	} else {
		bool b_even = vargs[1].representsEven();
		for(size_t i = 0; i < mstruct.size(); i++) {
			if(!b_even || !mstruct[i].representsReal(true)) mstruct[i].transformById(FUNCTION_ID_ABS);
			mstruct[i] ^= vargs[1];
		}
		if(mstruct.size() == 0) mstruct.clear(true);
		else if(mstruct.size() == 1) mstruct.setToChild(1, true);
		else mstruct.setType(STRUCT_ADDITION);
		mstruct ^= vargs[1];
		mstruct.last().inverse();
	}
	return 1;
}
bool contains_nonlinear_unit(const MathStructure &m) {
	if(m.isUnit()) return m.unit()->hasNonlinearRelationToBase();
	for(size_t i = 0; i < m.size(); i++) {
		if(contains_nonlinear_unit(m[i])) return true;
	}
	return false;
}
bool matrix_to_rref(MathStructure &m, const EvaluationOptions &eo2) {
	if(contains_nonlinear_unit(m)) return false;
	size_t rows = m.rows();
	size_t cols = m.columns();
	EvaluationOptions eo = eo2;
	eo.keep_zero_units = false;
	size_t cur_row = 0;
	for(size_t c = 0; c < cols; ) {
		bool b = false;
		for(size_t r = cur_row; r < rows; r++) {
			if(m[r][c].representsNonZero(true)) {
				if(r != cur_row) {
					MathStructure *mrow = &m[r];
					mrow->ref();
					m.delChild(r + 1);
					m.insertChild_nocopy(mrow, cur_row + 1);
				}
				for(r = 0; r < rows; r++) {
					if(r != cur_row) {
						if(m[r][c].representsNonZero(true)) {
							MathStructure mmul(m[r][c]);
							mmul.calculateDivide(m[cur_row][c], eo);
							mmul.calculateNegate(eo);
							for(size_t c2 = 0; c2 < cols; c2++) {
								if(c2 == c) {
									m[r][c2].clear(true);
								} else {
									MathStructure madd(m[cur_row][c2]);
									madd.calculateMultiply(mmul, eo);
									m[r][c2].calculateAdd(madd, eo);
								}
							}
						} else if(!m[r][c].representsZero(true)) {
							return false;
						}
					}
				}
				for(size_t c2 = 0; c2 < cols; c2++) {
					if(c2 != c) {
						m[cur_row][c2].calculateDivide(m[cur_row][c], eo);
					}
				}
				m[cur_row][c].set(1, 1, 0, true);
				cur_row++;
				b = true;
				break;
			} else if(!m[r][c].representsZero(true)) {
				return false;
			}
		}
		if(cur_row == rows) break;
		if(!b) c++;
	}
	return true;
}
RRefFunction::RRefFunction() : MathFunction("rref", 1) {
	setArgumentDefinition(1, new MatrixArgument());
}
int RRefFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	// echelon matrix
	MathStructure m(vargs[0]);
	if(!matrix_to_rref(m, eo)) return false;
	mstruct = m;
	return 1;
}
MatrixRankFunction::MatrixRankFunction() : MathFunction("rk", 1) {
	setArgumentDefinition(1, new MatrixArgument());
}
int MatrixRankFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	MathStructure m(vargs[0]);
	if(!matrix_to_rref(m, eo)) return false;
	size_t rows = m.rows();
	size_t cols = m.columns();
	Number nr;
	// count zero rows
	for(size_t r = 0; r < rows; r++) {
		bool b = false;
		for(size_t c = 0; c < cols; c++) {
			if(m[r][c].representsNonZero(true)) {
				b = true;
				break;
			} else if(!m[r][c].representsZero(true)) {
				return false;
			}
		}
		if(!b) break;
		nr++;
	}
	mstruct = nr;
	return 1;
}
EntrywiseFunction::EntrywiseFunction() : MathFunction("entrywise", 2) {
	VectorArgument *arg = new VectorArgument();
	arg->addArgument(new VectorArgument());
	arg->addArgument(new SymbolicArgument());
	arg->setReoccuringArguments(true);
	setArgumentDefinition(2, arg);
}
int EntrywiseFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	if(vargs[1].size() == 0) {
		mstruct = vargs[0];
		return 1;
	}
	bool b_matrix = vargs[1][0].isMatrix();
	for(size_t i3 = 2; i3 < vargs[1].size(); i3 += 2) {
		if(b_matrix) {
			if(!vargs[1][i3].isMatrix() || vargs[1][i3].columns() != vargs[1][0].columns() || vargs[1][i3].rows() != vargs[1][0].rows()) {
				CALCULATOR->error(true, _("%s() requires that all matrices/vectors have the same dimensions."), preferredDisplayName().name.c_str(), NULL);
				return 0;
			}
		} else if(vargs[1][i3].size() != vargs[1][0].size()) {
			CALCULATOR->error(true, _("%s() requires that all matrices/vectors have the same dimensions."), preferredDisplayName().name.c_str(), NULL);
			return 0;
		}
	}
	MathStructure mexpr(vargs[0]);
	EvaluationOptions eo2 = eo;
	eo2.calculate_functions = false;
	mexpr.eval(eo2);
	mstruct = vargs[1][0];
	for(size_t i = 0; i < mstruct.size(); i++) {
		if(b_matrix) {
			for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
				mstruct[i][i2] = mexpr;
				for(size_t i3 = 1; i3 < vargs[1].size(); i3 += 2) {
					mstruct[i][i2].replace(vargs[1][i3], vargs[1][i3 - 1][i][i2]);
				}
			}
		} else {
			mstruct[i] = mexpr;
			for(size_t i3 = 1; i3 < vargs[1].size(); i3 += 2) {
				mstruct[i].replace(vargs[1][i3], vargs[1][i3 - 1][i]);
			}
		}
	}
	return 1;
}

DotProductFunction::DotProductFunction() : MathFunction("dot", 2) {
	setArgumentDefinition(1, new VectorArgument(""));
	setArgumentDefinition(2, new VectorArgument(""));
}
int DotProductFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct = vargs[0];
	MathStructure m2(vargs[1]);
	if(mstruct.size() == m2.size()) {
		for(size_t i = 0; i < mstruct.size(); i++) {
			mstruct[i] *= m2[i];
		}
		if(mstruct.size() == 1) mstruct.setToChild(1);
		else mstruct.setType(STRUCT_ADDITION);
		return 1;
	}
	return 0;
}

EntrywiseMultiplicationFunction::EntrywiseMultiplicationFunction() : MathFunction("times", 1, -1) {}
int EntrywiseMultiplicationFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	if(vargs.size() == 0) {mstruct.clear(); return 1;}
	if(vargs.size() == 1) {
		mstruct = vargs[0];
		if(!mstruct.isVector()) {
			if(mstruct.representsScalar()) return 1;
			CALCULATOR->beginTemporaryStopMessages();
			mstruct.eval(eo);
			if(!mstruct.isVector()) {
				CALCULATOR->endTemporaryStopMessages();
				mstruct = vargs[0];
				return 1;
			}
			CALCULATOR->endTemporaryStopMessages(true);
		}
		mstruct.setType(STRUCT_FUNCTION);
		mstruct.setFunction(this);
		return 1;
	}
	if(vargs.size() > 2) {
		mstruct = vargs[0];
		bool b_scalar = mstruct.representsScalar();
		for(size_t index = 1; index < vargs.size(); index++) {
			if(CALCULATOR->aborted()) return 0;
			if(b_scalar || vargs[index].representsScalar()) {
				mstruct.multiply(vargs[index], true);
				if(b_scalar) b_scalar = mstruct.representsScalar();
			} else {
				mstruct.transform(STRUCT_VECTOR, vargs[index]);
				mstruct.transform(this);
			}
		}
		return 1;
	}
	mstruct = vargs[0];
	MathStructure m2(vargs[1]);
	bool b_eval = false;
	if(mstruct.representsScalar() || m2.representsScalar()) {
		mstruct *= m2;
		return 1;
	}
	if(!mstruct.isVector() || (!mstruct.isMatrix() && !mstruct.representsNonMatrix())) {
		mstruct.eval(eo);
		if(mstruct.representsScalar()) {mstruct *= m2; return 1;}
		b_eval = true;
	}
	if(!m2.isVector() || (!m2.isMatrix() && !m2.representsNonMatrix())) {
		m2.eval(eo);
		if(m2.representsScalar()) {mstruct *= m2; return 1;}
		b_eval = true;
	}
	if(mstruct.isVector() && m2.isVector()) {
		if(mstruct.isMatrix()) {
			if(m2.isMatrix()) {
				if(mstruct.size() == m2.size()) {
					if(mstruct[0].size() == m2[0].size()) {
						for(size_t i = 0; i < mstruct.size(); i++) {
							for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
								if(CALCULATOR->aborted()) return 0;
								mstruct[i][i2] *= m2[i][i2];
							}
						}
						return 1;
					} else if(mstruct[0].size() == 1) {
						for(size_t i = 0; i < m2.size(); i++) {
							for(size_t i2 = 1; i2 < m2[i].size(); i2++) {
								if(CALCULATOR->aborted()) return 0;
								mstruct[i].addChild(mstruct[i][0]);
							}
							for(size_t i2 = 0; i2 < m2[i].size(); i2++) {
								if(CALCULATOR->aborted()) return 0;
								mstruct[i][i2] *= m2[i][i2];
							}
						}
						return 1;
					} else if(m2[0].size() == 1) {
						for(size_t i = 0; i < mstruct.size(); i++) {
							for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
								if(CALCULATOR->aborted()) return 0;
								mstruct[i][i2] *= m2[i][0];
							}
						}
						return 1;
					}
				}
			} else if(mstruct.columns() == 1) {
				for(size_t i = 0; i < mstruct.size(); i++) {
					for(size_t i2 = 1; i2 < m2.size(); i2++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i].addChild(mstruct[i][0]);
					}
					for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i][i2] *= m2[i2];
					}
				}
				return 1;
			}
		} else {
			if(m2.isMatrix() && m2.columns() == 1) {
				mstruct.transform(STRUCT_VECTOR);
				for(size_t i = 1; i < m2.size(); i++) {
					if(CALCULATOR->aborted()) return 0;
					mstruct.addChild(mstruct[0]);
				}
				for(size_t i = 0; i < mstruct.size(); i++) {
					for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i][i2] *= m2[i][0];
					}
				}
				return 1;
			} else if(mstruct.size() == m2.size()) {
				for(size_t i = 0; i < mstruct.size(); i++) {
					if(CALCULATOR->aborted()) return 0;
					mstruct[i] *= m2[i];
				}
				return 1;
			}
		}
	}
	if(!b_eval) return 0;
	mstruct.transform(STRUCT_VECTOR, m2);
	return -1;
}
EntrywiseDivisionFunction::EntrywiseDivisionFunction() : MathFunction("rdivide", 2) {
}
int EntrywiseDivisionFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct = vargs[0];
	MathStructure m2(vargs[1]);
	bool b_eval = false;
	if(m2.representsScalar()) {
		mstruct /= m2;
		return 1;
	}
	if(!mstruct.isVector() || (!mstruct.isMatrix() && !mstruct.representsNonMatrix())) {
		mstruct.eval(eo);
		b_eval = true;
	}
	if(!m2.isVector() || (!m2.isMatrix() && !m2.representsNonMatrix())) {
		m2.eval(eo);
		if(m2.representsScalar()) {mstruct /= m2; return 1;}
		b_eval = true;
	}
	if(m2.isVector()) {
		if(mstruct.representsScalar()) {
			MathStructure m(mstruct);
			mstruct.clearVector();
			if(m2.isMatrix()) {
				mstruct.resizeVector(m2.size(), m_zero);
				if(mstruct.size() < m2.size()) return 0;
				for(size_t i = 0; i < m2.size(); i++) {
					mstruct[i].clearVector();
					mstruct[i].resizeVector(m2[i].size(), m);
					if(mstruct[i].size() < m2[i].size()) return 0;
					for(size_t i2 = 0; i2 < m2[i].size(); i2++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i][i2] = m;
						mstruct[i][i2] /= m2[i][i2];
					}
				}
				return 1;
			} else {
				mstruct.resizeVector(m2.size(), m);
				if(mstruct.size() < m2.size()) return 0;
				for(size_t i = 0; i < m2.size(); i++) {
					if(CALCULATOR->aborted()) return 0;
					mstruct[i] /= m2[i];
				}
				return 1;
			}
		}
		if(mstruct.isVector()) {
			if(mstruct.isMatrix()) {
				if(m2.isMatrix()) {
					if(m2.size() == mstruct.size()) {
						if(mstruct[0].size() == m2[0].size()) {
							for(size_t i = 0; i < mstruct.size(); i++) {
								for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i][i2] /= m2[i][i2];
								}
							}
							return 1;
						} else if(mstruct[0].size() == 1) {
							for(size_t i = 0; i < m2.size(); i++) {
								for(size_t i2 = 1; i2 < m2[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i].addChild(mstruct[i][0]);
								}
								for(size_t i2 = 0; i2 < m2[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i][i2] /= m2[i][i2];
								}
							}
							return 1;
						} else if(m2[0].size() == 1) {
							for(size_t i = 0; i < mstruct.size(); i++) {
								for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i][i2] /= m2[i][0];
								}
							}
							return 1;
						}
					}
				} else if(mstruct.columns() == 1) {
					for(size_t i = 0; i < mstruct.size(); i++) {
						for(size_t i2 = 1; i2 < m2.size(); i2++) {
							if(CALCULATOR->aborted()) return 0;
							mstruct[i].addChild(mstruct[i][0]);
						}
						for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
							if(CALCULATOR->aborted()) return 0;
							mstruct[i][i2] /= m2[i2];
						}
					}
					return 1;
				}
			} else {
				if(m2.isMatrix() && m2.columns() == 1) {
					mstruct.transform(STRUCT_VECTOR);
					for(size_t i = 1; i < m2.size(); i++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct.addChild(mstruct[0]);
					}
					for(size_t i = 0; i < mstruct.size(); i++) {
						for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
							if(CALCULATOR->aborted()) return 0;
							mstruct[i][i2] /= m2[i][0];
						}
					}
					return 1;
				} else if(mstruct.size() == m2.size()) {
					for(size_t i = 0; i < mstruct.size(); i++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i] /= m2[i];
					}
					return 1;
				}
			}
		}
	}
	if(!b_eval) return 0;
	mstruct.transform(STRUCT_VECTOR, m2);
	return -3;
}
EntrywisePowerFunction::EntrywisePowerFunction() : MathFunction("pow", 2) {
}
int EntrywisePowerFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	mstruct = vargs[0];
	MathStructure m2(vargs[1]);
	bool b_eval = false;
	if(!mstruct.representsScalar() && (!mstruct.isVector() || (!mstruct.isMatrix() && !mstruct.representsNonMatrix()))) {
		mstruct.eval(eo);
		b_eval = true;
	}
	if(!m2.representsScalar() && (!m2.isVector() || (!m2.isMatrix() && !m2.representsNonMatrix()))) {
		m2.eval(eo);
		b_eval = true;
	}
	if(m2.representsScalar()) {
		if(mstruct.representsScalar()) {
			mstruct ^= m2;
			return 1;
		} else if(mstruct.isVector()) {
			if(mstruct.isMatrix()) {
				for(size_t i = 0; i < mstruct.size(); i++) {
					for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i][i2] ^= m2;
					}
				}
				return 1;
			} else {
				for(size_t i = 0; i < mstruct.size(); i++) {
					if(CALCULATOR->aborted()) return 0;
					mstruct[i] ^= m2;
				}
				return 1;
			}
		}
	}
	if(m2.isVector()) {
		if(mstruct.representsScalar()) {
			MathStructure m(mstruct);
			mstruct.clearVector();
			if(m2.isMatrix()) {
				mstruct.resizeVector(m2.size(), m_zero);
				if(mstruct.size() < m2.size()) return 0;
				for(size_t i = 0; i < m2.size(); i++) {
					mstruct[i].clearVector();
					mstruct[i].resizeVector(m2[i].size(), m);
					if(mstruct[i].size() < m2[i].size()) return 0;
					for(size_t i2 = 0; i2 < m2[i].size(); i2++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i][i2] = m;
						mstruct[i][i2] ^= m2[i][i2];
					}
				}
				return 1;
			} else {
				mstruct.resizeVector(m2.size(), m);
				if(mstruct.size() < m2.size()) return 0;
				for(size_t i = 0; i < m2.size(); i++) {
					if(CALCULATOR->aborted()) return 0;
					mstruct[i] ^= m2[i];
				}
				return 1;
			}
		} else if(mstruct.isVector()) {
			if(mstruct.isMatrix()) {
				if(m2.isMatrix()) {
					if(m2.size() == mstruct.size()) {
						if(mstruct[0].size() == m2[0].size()) {
							for(size_t i = 0; i < mstruct.size(); i++) {
								for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i][i2] ^= m2[i][i2];
								}
							}
							return 1;
						} else if(mstruct[0].size() == 1) {
							for(size_t i = 0; i < m2.size(); i++) {
								for(size_t i2 = 1; i2 < m2[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i].addChild(mstruct[i][0]);
								}
								for(size_t i2 = 0; i2 < m2[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i][i2] ^= m2[i][i2];
								}
							}
							return 1;
						} else if(m2[0].size() == 1) {
							for(size_t i = 0; i < mstruct.size(); i++) {
								for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
									if(CALCULATOR->aborted()) return 0;
									mstruct[i][i2] ^= m2[i][0];
								}
							}
							return 1;
						}
					}
				} else if(mstruct.columns() == 1) {
					for(size_t i = 0; i < mstruct.size(); i++) {
						for(size_t i2 = 1; i2 < m2.size(); i2++) {
							if(CALCULATOR->aborted()) return 0;
							mstruct[i].addChild(mstruct[i][0]);
						}
						for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
							if(CALCULATOR->aborted()) return 0;
							mstruct[i][i2] ^= m2[i2];
						}
					}
					return 1;
				}
			} else {
				if(m2.isMatrix() && m2.columns() == 1) {
					mstruct.transform(STRUCT_VECTOR);
					for(size_t i = 1; i < m2.size(); i++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct.addChild(mstruct[0]);
					}
					for(size_t i = 0; i < mstruct.size(); i++) {
						for(size_t i2 = 0; i2 < mstruct[i].size(); i2++) {
							if(CALCULATOR->aborted()) return 0;
							mstruct[i][i2] ^= m2[i][0];
						}
					}
					return 1;
				} else if(mstruct.size() == m2.size()) {
					for(size_t i = 0; i < mstruct.size(); i++) {
						if(CALCULATOR->aborted()) return 0;
						mstruct[i] ^= m2[i];
					}
					return 1;
				}
			}
		}
	}
	if(!b_eval) return 0;
	mstruct.transform(STRUCT_VECTOR, m2);
	return -3;
}

bool process_replace(MathStructure &mprocess, const MathStructure &mstruct, const MathStructure &vargs, size_t index);
bool process_replace(MathStructure &mprocess, const MathStructure &mstruct, const MathStructure &vargs, size_t index) {
	if(mprocess == vargs[1]) {
		mprocess = mstruct[index];
		return true;
	}
	if(!vargs[3].isEmptySymbol() && mprocess == vargs[3]) {
		mprocess = (int) index + 1;
		return true;
	}
	if(!vargs[4].isEmptySymbol() && mprocess == vargs[4]) {
		mprocess = vargs[2];
		return true;
	}
	bool b = false;
	for(size_t i = 0; i < mprocess.size(); i++) {
		if(process_replace(mprocess[i], mstruct, vargs, index)) {
			mprocess.childUpdated(i + 1);
			b = true;
		}
	}
	return b;
}

ProcessFunction::ProcessFunction() : MathFunction("process", 3, 5) {
	setArgumentDefinition(2, new SymbolicArgument());
	setArgumentDefinition(3, new VectorArgument());
	setArgumentDefinition(4, new SymbolicArgument());
	setDefaultValue(4, "\"\"");
	setArgumentDefinition(5, new SymbolicArgument());
	setDefaultValue(5, "\"\"");
}
int ProcessFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {

	mstruct = vargs[2];
	MathStructure mprocess;
	for(size_t index = 0; index < mstruct.size(); index++) {
		mprocess = vargs[0];
		process_replace(mprocess, mstruct, vargs, index);
		mstruct[index] = mprocess;
	}
	return 1;

}

bool process_matrix_replace(MathStructure &mprocess, const MathStructure &mstruct, const MathStructure &vargs, size_t rindex, size_t cindex);
bool process_matrix_replace(MathStructure &mprocess, const MathStructure &mstruct, const MathStructure &vargs, size_t rindex, size_t cindex) {
	if(mprocess == vargs[1]) {
		mprocess = mstruct[rindex][cindex];
		return true;
	}
	if(!vargs[3].isEmptySymbol() && mprocess == vargs[3]) {
		mprocess = (int) rindex + 1;
		return true;
	}
	if(!vargs[4].isEmptySymbol() && mprocess == vargs[4]) {
		mprocess = (int) cindex + 1;
		return true;
	}
	if(!vargs[5].isEmptySymbol() && mprocess == vargs[5]) {
		mprocess = vargs[2];
		return true;
	}
	bool b = false;
	for(size_t i = 0; i < mprocess.size(); i++) {
		if(CALCULATOR->aborted()) break;
		if(process_matrix_replace(mprocess[i], mstruct, vargs, rindex, cindex)) {
			mprocess.childUpdated(i + 1);
			b = true;
		}
	}
	return b;
}

ProcessMatrixFunction::ProcessMatrixFunction() : MathFunction("processm", 3, 6) {
	setArgumentDefinition(2, new SymbolicArgument());
	setArgumentDefinition(3, new MatrixArgument());
	setArgumentDefinition(4, new SymbolicArgument());
	setDefaultValue(4, "\"\"");
	setArgumentDefinition(5, new SymbolicArgument());
	setDefaultValue(5, "\"\"");
	setArgumentDefinition(6, new SymbolicArgument());
	setDefaultValue(6, "\"\"");
}
int ProcessMatrixFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {

	mstruct = vargs[2];
	MathStructure mprocess;
	for(size_t rindex = 0; rindex < mstruct.size(); rindex++) {
		for(size_t cindex = 0; cindex < mstruct[rindex].size(); cindex++) {
			mprocess = vargs[0];
			process_matrix_replace(mprocess, mstruct, vargs, rindex, cindex);
			if(CALCULATOR->aborted()) return 0;
			mstruct[rindex][cindex] = mprocess;
		}
	}
	return 1;

}
GenerateVectorFunction::GenerateVectorFunction() : MathFunction("genvector", 3, 6) {
	setDefaultValue(4, "1");
	setArgumentDefinition(5, new SymbolicArgument());
	setDefaultValue(5, "undefined");
	IntegerArgument *iarg = new IntegerArgument("");
	iarg->setMax(&nr_one);
	iarg->setMin(&nr_minus_one);
	setArgumentDefinition(6, iarg);
	setDefaultValue(6, "-1");
}
int GenerateVectorFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	if(CALCULATOR->aborted()) return 0;
	bool b_step = vargs[5].number().isPositive();
	MathStructure msteps(vargs[3]);
	if(!b_step) {
		CALCULATOR->beginTemporaryStopMessages();
		msteps.eval(eo);
		if(!msteps.number().isInteger() && msteps.number().isInterval() && !msteps.number().isNonInteger() && msteps.number().precision(true) > PRECISION + 10) {
			Number nr;
			if(msteps.number().getCentralInteger(nr)) {
				msteps.set(nr, true);
			}
		}
		if(vargs[5].number().isNegative()) b_step = !msteps.isInteger() || msteps.number().isNegative() || msteps.number().isOne();
		CALCULATOR->endTemporaryStopMessages(!b_step);
	}
	if(b_step) {
		mstruct = vargs[0].generateVector(vargs[4], vargs[1], vargs[2], vargs[3], NULL, eo);
	} else {
		bool overflow = false;
		int steps = msteps.number().intValue(&overflow);
		if(!msteps.isInteger() || overflow || steps < 1) {
			CALCULATOR->error(true, _("The number of requested elements in generate vector function must be a positive integer."), NULL);
			return 0;
		}
		mstruct = vargs[0].generateVector(vargs[4], vargs[1], vargs[2], steps, NULL, eo);
	}
	if(CALCULATOR->aborted() || mstruct.size() == 0) return 0;
	return 1;
}
ColonFunction::ColonFunction() : MathFunction("colon", 2, 3) {
	Argument *arg = new Argument();
	arg->setHandleVector(true);
	setArgumentDefinition(1, arg);
	arg = new Argument();
	arg->setHandleVector(true);
	setArgumentDefinition(2, arg);
	arg = new Argument();
	arg->setHandleVector(true);
	setArgumentDefinition(3, arg);
	setDefaultValue(3, "undefined");
}
int ColonFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	if(CALCULATOR->aborted()) return 0;
	mstruct.set("x", true);
	if(vargs[2].isUndefined()) mstruct = mstruct.generateVector(mstruct, vargs[0], vargs[1], m_one, NULL, eo);
	else mstruct = mstruct.generateVector(mstruct, vargs[0], vargs[2], vargs[1], NULL, eo);
	if(CALCULATOR->aborted() || mstruct.size() == 0) return 0;
	return 1;
}
SelectFunction::SelectFunction() : MathFunction("select", 2, 4) {
	setArgumentDefinition(1, new VectorArgument());
	setArgumentDefinition(3, new SymbolicArgument());
	setDefaultValue(3, "undefined");
	setArgumentDefinition(4, new BooleanArgument());
	setDefaultValue(4, "0");
}
int SelectFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions &eo) {
	MathStructure mtest;
	mstruct = vargs[0];
	for(size_t i = 0; i < mstruct.size();) {
		if(CALCULATOR->aborted()) return 0;
		mtest = vargs[1];
		mtest.replace(vargs[2], mstruct[i], false, false, true);
		mtest.eval(eo);
		if(!mtest.isNumber() || mtest.number().getBoolean() < 0) {
			CALCULATOR->error(true, _("Comparison failed."), NULL);
			return 0;
		}
		if(mtest.number().getBoolean() == 0) {
			if(vargs[3].number().getBoolean() == 0) {
				mstruct.delChild(i + 1);
			} else {
				i++;
			}
		} else if(vargs[3].number().getBoolean() > 0) {
			mstruct.setToChild(i + 1);
			return 1;
		} else {
			i++;
		}
	}
	if(vargs[3].number().getBoolean() > 0) {
		CALCULATOR->error(true, _("No matching item found."), NULL);
		return 0;
	}
	return 1;
}
LoadFunction::LoadFunction() : MathFunction("load", 1, 3) {
	setArgumentDefinition(1, new FileArgument());
	Argument *arg = new IntegerArgument("", ARGUMENT_MIN_MAX_POSITIVE, true, true, INTEGER_TYPE_SINT);
	arg->setHandleVector(false);
	setArgumentDefinition(2, arg);
	setDefaultValue(2, "1");
	setArgumentDefinition(3, new TextArgument());
	setDefaultValue(3, ",");
}
int LoadFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	string delim = vargs[2].symbol();
	if(delim == "tab") {
		delim = "\t";
	}
	if(!CALCULATOR->importCSV(mstruct, vargs[0].symbol().c_str(), vargs[1].number().intValue(), delim)) {
		CALCULATOR->error(true, "Failed to load %s.", vargs[0].symbol().c_str(), NULL);
		return 0;
	}
	return 1;
}
ExportFunction::ExportFunction() : MathFunction("export", 2, 3) {
	setArgumentDefinition(1, new MatrixArgument());
	setArgumentDefinition(2, new FileArgument());
	setArgumentDefinition(3, new TextArgument());
	setDefaultValue(3, ",");
}
int ExportFunction::calculate(MathStructure&, const MathStructure &vargs, const EvaluationOptions&) {
	string delim = vargs[2].symbol();
	if(delim == "tab") {
		delim = "\t";
	}
	if(!CALCULATOR->exportCSV(vargs[0], vargs[1].symbol().c_str(), delim)) {
		CALCULATOR->error(true, "Failed to export to %s.", vargs[1].symbol().c_str(), NULL);
		return 0;
	}
	return 1;
}
KroneckerProductFunction::KroneckerProductFunction() : MathFunction("kron", 2) {
	setArgumentDefinition(1, new MatrixArgument());
	setArgumentDefinition(2, new MatrixArgument());
}
int KroneckerProductFunction::calculate(MathStructure &mstruct, const MathStructure &vargs, const EvaluationOptions&) {
	size_t r1 = vargs[0].rows(), r2 = vargs[1].rows(), c1 = vargs[0].columns(), c2 = vargs[1].columns();
	if(r2 > ((size_t) -1) / r1 || c2 > ((size_t) -1) / c1) return 0;
	mstruct.clearMatrix();
	mstruct.resizeMatrix(r1 * r2, c1 * c2, m_zero);
	size_t r, c;
	for(size_t ra = 0; ra < r1; ra++) {
		for(size_t ca = 0; ca < c1; ca++) {
			for(size_t rb = 0; rb < r2; rb++) {
				for(size_t cb = 0; cb < c2; cb++) {
					if(CALCULATOR->aborted()) return 0;
					r = ra * r2 + rb;
					c = ca * c2 + cb;
					mstruct[r][c] = vargs[0][ra][ca];
					mstruct[r][c] *= vargs[1][rb][cb];
				}
			}
		}
	}
	return 1;
}
