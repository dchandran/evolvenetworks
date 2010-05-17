#ifndef JSIMSBWOPTIMGUI_H
#define JSIMSBWOPTIMGUI_H

#include <QApplication>
#include <QString>
#include <QStringList>
#include <QWidget>
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QLabel>
#include <QToolBar>
#include <QMainWindow>
#include <QTableWidget>
#include <QCheckBox>
#include <QDockWidget>
#include <QDir>
#include <QMessageBox>
#include <QComboBox>
#include <QThread>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QGridLayout>
#include <QTextEdit>
#include <vector>
#include <iostream>
#include "SBW/sbw.h"
#include "plotWidget.h"
#include "randomize_widget.h"

class JSimRun : public QThread
{
	Q_OBJECT

signals:

	void errorChanged(double);
	void resultsChanged(std::vector<double>);

public:

	std::string sbml;
	std::string csvfile;
	std::vector<std::string> paramNames;
	std::vector<std::string> varNames;
	std::vector<double> paramInitialValues;
	std::vector<std::string> optimParamNames;
	std::vector<double> optimParamValues;

	JSimRun(QObject * parent = 0): QThread(parent) {}
	Module jsimModule;
	void run();
};

class JSimOptimGUI : public QMainWindow
{
	Q_OBJECT

public:

	JSimOptimGUI(QWidget * parent = 0);
	bool connectToRoadRunner();
	bool connectToJsim();
	void disconnectJsim();

private:

	Module roadRunnerModule;
	Module jsimOptimizerModule;
	void initializeWidgets();

private:

	std::string sbmlFile;
	std::string dataFile;
	bool roadRunnerLoaded;
	bool jsimLoaded;
	double timeEnd, stepSize;
	int numDataColumns;

	DataPlot * errorPlot;

	DataPlot * initialSimPlot;
	DataPlot * currentSimPlot;

	DataPlot * residualPlot;
	DataPlot * dataFitPlot;

	QWidget * errorPlotWindow; //error vs iterations (1 plot)
	QWidget * simPlotWindow;  //initial simulation, current simulation (2 plots)
	QWidget * fittedPlotWindow; //data + current simulated data, residual (2 plots)

	int runCount, iterCount;

	QLabel * sbmlFileLabel;
	QLabel * dataFileLabel;

	QTableWidget * paramsTable;
	QTableWidget * advancedTable;
	QList<QCheckBox*> checkBoxes;
	QComboBox * algorithmsListBox;

	std::vector<std::string> paramNames;
	std::vector<double> paramValues;

	std::vector<std::string> optimizerParamNames;
	std::vector<double> optimizerParamValues;

	std::vector<std::string> varNames;
	std::vector<double> initialValues;

	QSize plotWindowSize;

	void fillParamsTable(const std::vector<double>& , const QStringList& );
	DataTable<double> simulate();
	void createGradedDataPlot(DataTable<qreal>&);

	RandomizerWidget * randomizer;
	QTextEdit * covMat;
	QLabel * errorLabel;

signals:
	
	void stop();

private slots:

	void randomizeInitialValues();
	void paramsTableChanged(int,int);
	void advancedTableChanged(int,int);
	void loadSBML();
	void loadData();
	void run();
	void algorithmChanged(const QString&);
	void updateErrorPlot(double);
	void updateOptimResults(std::vector<double>);
	void setErrorPlotVisible(bool);
	void setSimPlotVisible(bool);
	void setFitPlotVisible(bool);

protected slots:
	void keyPressEvent ( QKeyEvent * event );
	void mousePressEvent ( QMouseEvent * event );
};

#endif
