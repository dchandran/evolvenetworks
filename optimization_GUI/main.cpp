#include <iostream>
#include <vector>
#include "sbwclient.h"

using namespace std;

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

	qRegisterMetaType< std::vector<double> >("std::vector<double>");
	qRegisterMetaType< std::vector< std::vector<double> > >("std::vector< std::vector<double> >");
	//qRegisterMetaType< const std::vector<double> >("const std::vector<double>");

	JSimOptimGUI gui;
	gui.show();

	int res = app.exec();

	gui.disconnectJsim();

	return res;
}
