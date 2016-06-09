#ifndef TORSIONRULE_H_
#define TORSIONRULE_H_

#include <iostream>
#include <vector>

using namespace std;
class TorsionRule
{
public:
	TorsionRule(string piName, string piPattern, vector<unsigned int> piPoints, vector<double> piAngles);
	virtual ~TorsionRule();
	// Methods
	string getName();	
	string getPattern();	
	vector<unsigned int> getPoints();	
	vector<double> getAngles();	
private:
	// Attributes
	string name;	
	string pattern;
	vector<unsigned int> points;
	vector<double> angles;
};

#endif /*TORSIONRULE_H_*/
