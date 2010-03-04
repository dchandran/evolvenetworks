/****************************************************************************

 Copyright (C) 2008 Deepak Chandran
 Contact: Deepak Chandran (dchandran1@gmail.com)
 See COPYWRITE.TXT

 This is the main application file for Tinkercell. It constructs a MainWindow
 and loads a list of default plugins.

****************************************************************************/
#include <QApplication>
#include "DesignFFNWidget.h"

int main(int argc, char *argv[])
{
    QApplication::setColorSpec (QApplication::ManyColor);
    QApplication app(argc, argv);
    
    FeedforwardNetworkDesigner::Table widget;    
	widget.show();
	
	int output = app.exec();
	
	app.closeAllWindows();

    return output;
}

