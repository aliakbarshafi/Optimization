// -------------------------------------------------------------- -*- C++ -*-
// File: UFL_ALTERNATE.cpp
// Version 12.1.0  
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5724-Y48
// (c) Copyright IBM Corporation 2000, 2009. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
// 
// The problem is to determine the optimal flow from facilities to clients for a single product
// The Code has been tested for the sample input taken from UFL.dat 
// The parameters for model are:
// demand of clients (d),
// unit distribution cost from facilities to clients (c) and 
// Fixed cost of producing at facilities (f).
//
//
// The decision variables of the model are:
//   x[p][q]  : portion of demand of client q satisfied by facility p
//   y[p] = 1 if production occurs at facility p, p belongs to P and y[p]=0 otherwise.

// The objective function 
// minimize  sum(over all facilities) sum (over all clients)(c[p][q] * x[p][q])
//
// The constraints are
// sum(over p) x[p][q] = d[q]  for all q belong to Q
// x[p][q] <= d[q])*y[p]  for all q belong to Q and p belong to P
// 
//  The bounds on the variables are:
//			x[p][q] >= 0
//			y[p] =  {0,1}
//

#include <ilcplex/ilocplex.h>;
#include <fstream>;
#include <cstring>;
#include <fstream>;
#include <string>;
#include <vector>;
#include <functional>;
#include <iterator>;

ILOSTLBEGIN

//Objects initiated for 2 dimensional Integer Array for creating ARC variables.

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloBoolVarArray IloBoolVarArray;


/* Function to iterate input file */
bool iterateFile(std::string fileName, std::function<void(const std::string &)> callback)
{
	// Open the File
	std::ifstream in(fileName.c_str());

	// Check if object is valid
	if (!in)
	{
		std::cerr << "Cannot open the File : " << fileName << std::endl;
		return false;
	}

	std::string str;
	// Read the next line from File until it reaches the end.
	while (std::getline(in, str))
	{
		// Skip the blank lines
		if (str != "")
			callback(str);
	}
	in.close();
	return true;
}

/* Function: parseFile
Input: vector Array of Unparsed File
Output: File in Correct Input format required for CPLEX
*/

ifstream parseFile(vector<string> &vecOfStr, int argc, char **argv) {
	// Parsing and getting input file in correct format
	int length = vecOfStr.size();
	std::ofstream tempFile;
	tempFile.open("temp.txt", fstream::app);
	tempFile << vecOfStr[0];

	for (int i = 1; i < length - 2; i++) {
		std::size_t found = vecOfStr[i].find_last_of(",");
		tempFile << '\n' << '[' << vecOfStr[i].substr(0, found) << ']' << ",";
	}

	tempFile << '\n' << '[' << vecOfStr[length - 2] << ']';
	tempFile << '\n' << vecOfStr[length - 1];
	tempFile.close();

	const char* filename1 = "temp.txt";
	if (argc >= 2) filename1 = argv[1];
	ifstream file1(filename1);

	if (!file1) {
		cerr << "No such file: " << filename1 << endl;
		throw(-1);
	}
	return file1;
}


int main(int argc, char **argv)
{
	IloEnv env;
	fstream outputFile;
	try {

		// Reading the input file

		const char* filename = "UFL.dat";
		if (argc >= 2) filename = argv[1];
		ifstream file(filename);
		if (!file) {
			cerr << "No such file: " << filename << endl;
			throw(-1);
		}

		std::vector<std::string> vecOfStr;
		bool res = iterateFile("UFL.dat", [&](const std::string & str) {
			vecOfStr.push_back(str);
		});

		// Attachng the Parameters of model with environment
		IloNumArray2 c(env);
		IloNumArray f(env);
		IloNumArray d(env);

		ifstream file1 = parseFile(vecOfStr, argc, &*argv);
		file1 >> d >> c >> f;
		file1.close();
		std::remove("temp.txt");

		// Calculating number of Facilities and Clients.
		IloInt n_p = f.getSize();
		IloInt n_q = d.getSize();
	
		/*Defining indices for facility and client
		p: facility
		q: client
		*/
		IloInt p, q;

		/* Creaing two models for MIP and LP relaxation
		mod: MIP problem
		modLP: LP Relaxation
		*/
		IloModel mod(env);
		IloModel modLP(env);

		// Attach Variables to model for each combination of facility and client
		IloNumVarArray2 x(env);
		IloBoolVarArray y(env);

		for (p = 0; p < n_p; p++) {
			x.add(IloNumVarArray(env, n_q, 0.0, IloInfinity));
		}
		
		y.add(IloBoolVarArray(env, n_p));

		// Attach Objective Function Cost to models: 
		IloExpr Cost(env);
		for (p = 0; p < n_p; p++) {
			for (q = 0; q < n_q; q++) {
				Cost += c[p][q] * x[p][q];
			}
		}

		for (p = 0; p < n_p; p++) {
			Cost += f[p] * y[p];
		}
		mod.add(IloMinimize(env, Cost));
		modLP.add(IloMinimize(env, Cost));
		Cost.end();
		
		// Adding constraints to model

		// sum(over p) x[p][q] = d[q]  for all q belong to Q	
		for (q = 0; q < n_q; q++) {
			IloExpr demandConstraint(env);
			for (p = 0; p < n_p; p++) {
				demandConstraint += x[p][q];
			}
			mod.add(demandConstraint == d[q]);
			modLP.add(demandConstraint == d[q]);
			demandConstraint.end();
		}

		// x[p][q] <= d[q]*y[p]  for all q belong to Q and p belong to P
		for (p = 0; p < n_p; p++) {
			for (q = 0; q < n_q; q++) {
				mod.add(x[p][q] <= d[q] * y[p]);
				modLP.add(x[p][q] <= d[q] * y[p]);

			}
		}
		
		// Relaxing the binary constraint on y for LP Relaxation and attaching a cplex object to solve  MIP and LP relaxation
		modLP.add(IloConversion(env, y, ILOFLOAT));

		IloCplex cplex(mod);
		IloCplex cplexLP(modLP);

		// Setting default cuts OFF for cplex and cplexLP */
		cplex.setParam(IloCplex::EachCutLim, 0);
		cplexLP.setParam(IloCplex::EachCutLim, 0);

		env.out() << endl << "*****************************    Section E   *****************************" << endl << endl;
		
		cplex.setParam(IloCplex::ClockType, 1);
		cplexLP.setParam(IloCplex::ClockType, 1);
		env.out()<<cplex.getNintVars();
		cplex.solve();
		cplexLP.solve();
				
		env.out() << endl << "Optimal Objective Value of MIP instance: " << cplex.getObjValue() << endl;
		env.out() << endl << "No of Integer Variables: " << n_p<< endl;
		env.out() << endl << "No of Continuous Variables: " << n_p*n_q << endl;
		env.out() << endl << "No of Constraints: " << cplex.getNrows() << endl;
		env.out() << endl << "Run time to solve the LP relaxation: " << cplexLP.getTime() << endl;
		env.out() << endl << "Optimal objective function value for this LP relaxation: " << cplexLP.getObjValue() << endl;
		env.out() << endl << "Run Time to solve MIP: " << cplex.getTime() << endl;
		env.out() << endl << "No of Nodes: " << cplex.getNnodes() << endl;
		env.out() << endl << "Percentage Gap pf MIP and LP Solutions: " << 100 * (cplex.getObjValue() - cplexLP.getObjValue()) / cplex.getObjValue() << endl;
		env.out() << endl << "No of Cuts: " << cplex.getNcuts(IloCplex::CutMir) +
											   cplex.getNcuts(IloCplex::CutImplBd) +
											   cplex.getNcuts(IloCplex::CutFlowCover) +
											   cplex.getNcuts(IloCplex::CutFlowPath) << endl;
	}

	catch (IloException& ex) {
      cerr << "Error: " << ex << endl;
	}
	catch (...) {
      cerr << "Error: Unknown exception caught!" << endl;
	}

	env.end();
	int BP1; 
	
	return 0;
} // END main

