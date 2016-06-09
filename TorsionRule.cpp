#include "TorsionRule.h"

TorsionRule::TorsionRule(string piName, string piPattern, vector<unsigned int> piPoints, vector<double> piAngles)
{
	this->name = piName;
	this->pattern = piPattern;	
	this->points = piPoints;
	this->angles = piAngles;
}

TorsionRule::~TorsionRule()
{
}

// Methods

string TorsionRule::getName() {
	return this->name;
}

string TorsionRule::getPattern() {
	return this->pattern;
}

vector<unsigned int> TorsionRule::getPoints() {
	return this->points;
}

vector<double> TorsionRule::getAngles() {
	return this->angles;
}
