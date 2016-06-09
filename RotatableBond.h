#ifndef ROTATABLEBOND_H_
#define ROTATABLEBOND_H_

//! Includes & namespaces
#include "RotatableBond.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace OpenBabel;
using namespace std;

class RotatableBond
{
public:
	RotatableBond(string piType, vector<OBAtom*> piAtoms, vector<double> piAngles);
	virtual ~RotatableBond();
	// Methods
	string getType();
	vector<OBAtom*> getAtoms();
	vector<double> getAngles();
	void addAngle(double piAngle);
	void removeAngle(unsigned int piPosition);
	void setType(string piType);	
	void setAngles(vector<double> piAngles);	
	void updateAtom(unsigned int piPosition, OBAtom* piAtom);	
private:
	string type;
	vector<OBAtom*> atoms;
	vector<double> angles;
};

#endif /*ROTATABLEBOND_H_*/
