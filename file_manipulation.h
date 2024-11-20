#ifndef File_h
#define File_h

#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

template<typename T>
void Lettura_file(
	vector<T> & data_x,
	const string & fileName
)
{
	ifstream input_file ;
	input_file.open (fileName.c_str(), ios::in) ;
	
	T valueX;
	
	while (true){
		input_file>>valueX;
		
		data_x.push_back(valueX);

		if (input_file.eof () == true) break ;
	}
	
	input_file.close ();
};

template<typename T1, typename T2>
void Lettura_file(
	vector<T1> & data_x,
	vector<T2> & data_y,
	const string & fileName
)
{
	ifstream input_file ;
	input_file.open (fileName.c_str(), ios::in) ;
	
	T1 valueX;
	T2 valueY;
	
	while (true){
		input_file>>valueX;
		input_file>>valueY;
		
		data_x.push_back(valueX);
		data_y.push_back(valueY);

		if (input_file.eof () == true) break ;
	}
	
	input_file.close ();
};

template<typename T1, typename T2, typename T3>
void Lettura_file(
	vector<T1> & data_x,
	vector<T2> & data_y,
	vector<T3> & data_z,
	const string & fileName
)
{
	ifstream input_file ;
	input_file.open (fileName.c_str(), ios::in) ;
	
	T1 valueX;
	T2 valueY;
	T3 valueZ;
	
	while (true){
		input_file>>valueX;
		input_file>>valueY;
		input_file>>valueZ;
		
		data_x.push_back(valueX);
		data_y.push_back(valueY);
		data_z.push_back(valueZ);

		if (input_file.eof () == true) break ;
	}
	
	input_file.close ();
};

template<typename T>
void Scrittura_file(
	vector<T> & data_x,
	const string & fileName
)
{
	ofstream out_file;
	out_file.open (fileName.c_str(), ios::in);
	
	for(int i=0; i<data_x.size(); i++)	out_file<<data_x.at(i)<<"\n";
	
	out_file.close();
};

#endif
