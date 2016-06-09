#include "RotatableBond.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

RotatableBond::RotatableBond(string piType, vector<OBAtom*> piAtoms, vector<double> piAngles)
{
	this->type = piType;
	this->atoms = piAtoms;
	this->angles = piAngles;
}

RotatableBond::~RotatableBond()
{
}

// Methods

void RotatableBond::addAngle(double piAngle) {
	this->angles.push_back(piAngle);
}

void RotatableBond::removeAngle(unsigned int piPosition) {
	vector<double>::iterator itVector = this->angles.begin();
	for (unsigned int i = 0; i < angles.size(); i++) {
		if (i == piPosition) {
			this->angles.erase(itVector);
			break;
		} else {
			itVector++;
		}
	}
}

string RotatableBond::getType() {
	return type;
}

vector<OBAtom*> RotatableBond::getAtoms() {
	return atoms;
}

vector<double> RotatableBond::getAngles() {
	return angles;
}

void RotatableBond::setType(string piType) {
	this->type = piType;
}

void RotatableBond::setAngles(vector<double> piAngles) {
	this->angles = piAngles;
}

void RotatableBond::updateAtom(unsigned int piPosition, OBAtom* piAtom) {
	atoms[piPosition] = piAtom;
}
