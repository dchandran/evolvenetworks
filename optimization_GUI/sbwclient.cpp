#include <cstdlib>
#include <QFile>
#include <QToolButton>
#include <QDockWidget>
#include <QFileDialog>
//#include "JsimHeader.h"
#include "sbwclient.h"

using namespace std;
using namespace SystemsBiologyWorkbench;

void JSimRun::run()
{
	//JSim Optimization service
	try
	{
		Service optimService = jsimModule.findServiceByName("Optimization");
		ServiceMethod loadSBML = optimService.getMethod("void loadSBML(string)");
		ServiceMethod run = optimService.getMethod("void run()");
		ServiceMethod isRunning = optimService.getMethod("boolean isRunning()");
		ServiceMethod getResults = optimService.getMethod("double[] getResults()");
		ServiceMethod getError = optimService.getMethod("double getError()");
		ServiceMethod setTargetParameters = optimService.getMethod("void setTargetParameters(string[])");
		ServiceMethod setTargetParameterValues = optimService.getMethod("void setInitialParameterValues(double[])");
		ServiceMethod setTargetData = optimService.getMethod("void setTargetData(string)");
		ServiceMethod setTargetVariables = optimService.getMethod("void setTargetVariables(string[])");
		ServiceMethod setOptimizerParams = optimService.getMethod("void setOptimizationParameter(string, double)");

		loadSBML.call(DataBlockWriter() << sbml);
		setTargetParameters.call(DataBlockWriter() << paramNames);
		setTargetParameterValues.call(DataBlockWriter() << paramInitialValues);
		setTargetVariables.call(DataBlockWriter() << varNames);
		setTargetData.call(DataBlockWriter() << csvfile);

		for (int i=0; i < optimParamNames.size() && i < optimParamValues.size();++i)
		{
			setOptimizerParams.call(DataBlockWriter() << optimParamNames[i] << optimParamValues[i]);
		}

		DataBlockWriter emptyArgs;
		vector<double> results;
		double error;
		bool running = true;
		double oldError = 0.0;

		run.call(emptyArgs);
		while (running)
		{
			isRunning.call(emptyArgs) >> running;
			getResults.call(emptyArgs) >> results;
			getError.call(emptyArgs) >> error;

			if (oldError != error)
			{
				emit errorChanged(error);
				emit resultsChanged(results);
				oldError = error;
			}
		}

		getResults.call(emptyArgs) >> results;
		emit resultsChanged(results);
	}
	catch (SBWException * e)
	{
		std::cout << e->getMessage() << std::endl;
	}
}

JSimOptimGUI::JSimOptimGUI(QWidget * parent) : QMainWindow(parent)
{
	initializeWidgets();
	connectToJsim();
	runCount = 0;
	iterCount = 0;
	timeEnd = 100;
	stepSize = 0.1;
	numDataColumns = 0;
}

void JSimOptimGUI::randomizeInitialValues()
{
	int rows = paramsTable->rowCount();
	 
	for (int i=0; i < rows && i < initialValues.size(); ++i)
	{
		//double r = ((double)rand())/(double)RAND_MAX;
		//initialValues[i] = 2.0 * r * initialValues[i];
		initialValues[i] = randomizer->getRandomValue(i,initialValues[i]);
		paramsTable->setItem(i,1, new QTableWidgetItem(QString::number(initialValues[i])));
	}
}

void JSimOptimGUI::fillParamsTable(const vector<double>& vec, const QStringList& headers)
{
	paramsTable->setColumnCount(2);
	paramsTable->setRowCount(headers.size());
	paramsTable->setHorizontalHeaderLabels(QStringList() << "optimize?" << "initial value");
	paramsTable->setVerticalHeaderLabels(headers);
	paramNames.clear();
	
	for (int i=0; i < headers.size() && i < vec.size(); ++i)
	{
		QCheckBox * checkBox = new QCheckBox(paramsTable);
		
		checkBoxes << checkBox;
		paramNames.push_back(std::string(headers[i].toAscii().data()));

		paramsTable->setCellWidget(i,0,checkBox);
		paramsTable->setItem(i,1, new QTableWidgetItem(QString::number(vec[i])));
	}
	initialValues = vec;
}


DataTable<double> JSimOptimGUI::simulate()
{
	DataTable<double> dataTable;

	if (!roadRunnerLoaded) return dataTable;

	try
	{
		Service sim = roadRunnerModule.findServiceByName("sim");
		ServiceMethod simulate = sim.getMethod("double[][] simulate()");
		ServiceMethod setTimeStart = sim.getMethod("void setTimeStart(double)");
		ServiceMethod setTimeEnd = sim.getMethod("void setTimeEnd(double)");
		ServiceMethod setNumPoints = sim.getMethod("void setNumPoints(int)");
		ServiceMethod reset = sim.getMethod("void reset()");
		ServiceMethod getNames = sim.getMethod("{} getFloatingSpeciesNames()");

		DataBlockWriter emptyArgs;
		reset.call(emptyArgs);
		setTimeStart.call(DataBlockWriter() << 0);
		setTimeEnd.call(DataBlockWriter() << timeEnd);
		setNumPoints.call(DataBlockWriter() << (int)(timeEnd/stepSize));

		DataBlockReader dataReader = simulate.call(emptyArgs);
		double** rawData; 
		int numRows=0, numCols=0, i=0;
		dataReader.get(numRows, numCols,rawData);
	
		dataTable.resize(numRows,numCols);
		for (i=0; i < numRows; ++i)
			for (int j=0; j < numCols; ++j)
			{
				dataTable.value(i,j) = rawData[i][j];
			}
		DataBlockReader::free2DArray(numRows, rawData);

		DataBlockReader namesReader;
		getNames.call(emptyArgs) >> namesReader;
		
		i = 1;
		while (namesReader.getNextType() != SystemsBiologyWorkbench::TerminateType)
		{
		   string s;
		   namesReader >> s; 
		   dataTable.colName(i) = QString(s.c_str());
		   ++i;
		}

		dataTable.colName(0) = "time";
	}
	catch(SBWException *e)
	{
		std::cout << e->getMessage() << std::endl;
	}
	return dataTable;
}

bool JSimOptimGUI::connectToRoadRunner()
{
	roadRunnerLoaded = false;
	try
	{
		SBW::connect();
		roadRunnerModule = SBW::getModuleInstance("edu.kgi.roadRunner");
		
		Service sim = roadRunnerModule.findServiceByName("sim");
		ServiceMethod loadSBML = sim.getMethod("void loadSBMLFromFile(string)");
		ServiceMethod getParamNames = sim.getMethod("{} getGlobalParameterNames()");
		ServiceMethod getParamValues = sim.getMethod("double[] getGlobalParameterValues()");
		ServiceMethod getNumParams = sim.getMethod("int getNumberOfGlobalParameters()");
		
		DataBlockWriter emptyArgs;
		QStringList paramNames;
		vector<double> paramValues;
		int numParams;
		
		loadSBML.call(DataBlockWriter() << sbmlFile);
		getNumParams.call(emptyArgs) >> numParams;
		getParamValues.call(emptyArgs) >> paramValues;

		DataBlockReader paramReader;
		getParamNames.call(emptyArgs) >> paramReader;
		
		while (paramReader.getNextType() != SystemsBiologyWorkbench::TerminateType)
		{
			string s;
			paramReader >> s;
			paramNames << QString(s.c_str());
		}

		fillParamsTable(paramValues,paramNames);

		QList<double> paramValuesList;
		for (int i=0; i < paramValues.size(); ++i)
			paramValuesList << paramValues[i];

		randomizer->initialize(paramNames,paramValuesList);
		roadRunnerLoaded = true;

		DataTable<qreal> sim_results = simulate();
		initialSimPlot->plot(sim_results,0,"Initial Simulation");
		if (!dataFile.empty())
		{
			createGradedDataPlot(sim_results);
		}
	}
	catch(SBWException *e)
	{
		std::cout << e->getMessage() << std::endl;
	}
	return roadRunnerLoaded;
}

void JSimOptimGUI::algorithmChanged(const QString& name)
{
	if (!jsimLoaded) return;
	try
	{
		Service optimService = jsimOptimizerModule.findServiceByName("Optimization");
		ServiceMethod setAlg = optimService.getMethod("void setAlgorithm(string)");
		ServiceMethod isPointWtUsed = optimService.getMethod("boolean isPointWeightUsed()");
		ServiceMethod isUpBoundUsed = optimService.getMethod("boolean isUpperBoundUsed()");
		ServiceMethod isLowBoundUsed = optimService.getMethod("boolean isLowerBoundUsed()");
		ServiceMethod getParamNames = optimService.getMethod("string[] getOptimizationParameterNames()");
		ServiceMethod getParamValues = optimService.getMethod("string[] getOptimizationParameterValues()");

		setAlg.call(DataBlockWriter() << string(name.toAscii().data()));

		bool showHigh = false, showLow = false;
		vector<string> paramNames;
		vector<double> paramValues;

		isUpBoundUsed.call(DataBlockWriter()) >> showHigh;
		isLowBoundUsed.call(DataBlockWriter()) >> showLow;
		getParamNames.call(DataBlockWriter()) >> paramNames;
		getParamValues.call(DataBlockWriter()) >> paramValues;

		QStringList list;
		for (int i=0; i < paramNames.size(); ++i)
			list << QString(paramNames[i].c_str());

		advancedTable->setColumnCount(1);
		advancedTable->setRowCount(paramNames.size());
		advancedTable->setVerticalHeaderLabels(list);
		for (int i=0; i < list.size() && i < paramValues.size(); ++i)
		{
			advancedTable->setItem(i,0,new QTableWidgetItem(QString::number(paramValues[i])));
		}
		advancedTable->setHorizontalHeaderLabels(QStringList() << "value");

		if (showLow && showHigh)
		{
			paramsTable->setColumnCount(4);
			paramsTable->setHorizontalHeaderLabels(QStringList() << "optimize?" << "initial values" << "lower bound" << "upper bound");
		}
		else
		if (showLow || showHigh)
		{
			paramsTable->setColumnCount(3);
			if (showLow)
				paramsTable->setHorizontalHeaderLabels(QStringList() << "optimize?" << "initial values" << "lower bound");
			else
				paramsTable->setHorizontalHeaderLabels(QStringList() << "optimize?" << "initial values" << "upper bound");
		}
		else
		{
			paramsTable->setColumnCount(2);
			paramsTable->setHorizontalHeaderLabels(QStringList() << "optimize?" << "initial values");
		}
	}
	catch(SBWException * e)
	{
		std::cout << e->getMessage() << std::endl;
	}	
}

bool JSimOptimGUI::connectToJsim()
{
	jsimLoaded = false;
	try
	{
		SBW::connect();
		jsimOptimizerModule = SBW::getModuleInstance("JSimSBWOptimServer");
		Service optimService = jsimOptimizerModule.findServiceByName("Optimization");
		ServiceMethod jsimConnect = optimService.getMethod("void connect()");
		ServiceMethod allAlgorithms = optimService.getMethod("string[] getAlgorithmList()");

		jsimConnect.call(DataBlockWriter());

		vector<string> names;
		allAlgorithms.call(DataBlockWriter()) >> names;

		QStringList texts;
		for (int i=0; i < names.size(); ++i)
			texts << QString(names[i].c_str());

		while (algorithmsListBox->count() > 0)
			algorithmsListBox->removeItem(algorithmsListBox->count()-1);

		jsimLoaded = true;

		if (texts.size() > 0)
		{
			algorithmsListBox->addItems(texts);
			connect(algorithmsListBox,SIGNAL(currentIndexChanged ( const QString & )),
					this,SLOT(algorithmChanged( const QString & )));
			algorithmsListBox->setCurrentIndex(0);
		}
	}

	catch(SBWException *e)
	{
		std::cout << e->getMessage() << std::endl;
		
	}
	return jsimLoaded;
}

void JSimOptimGUI::initializeWidgets()
{
	QMenu * fileMenu = new QMenu("&File",this);
	QMenu * runMenu = new QMenu("&Run",this);
	QMenu * viewMenu = new QMenu("&View",this);
	QMenu * helpMenu = new QMenu("&Help",this);

	QAction * loadModelAction = fileMenu->addAction("&Load model",this,SLOT(loadSBML()));
	QAction * loadCSVAction = fileMenu->addAction("&Open data",this,SLOT(loadData()));
	QAction * exitAction = fileMenu->addAction("&Exit",this,SLOT(close()));

	QAction * runAction = runMenu->addAction("&Optimize",this,SLOT(run()));
	QAction * stopAction = runMenu->addAction("&Stop",this,SIGNAL(stop()));
	
	QToolBar * toolBar = new QToolBar(this);
	toolBar->addAction(loadModelAction);
	toolBar->addAction(loadCSVAction);
	toolBar->addAction(runAction);
	toolBar->addAction(stopAction);

	menuBar()->addMenu(fileMenu);
	menuBar()->addMenu(runMenu);
	menuBar()->addMenu(viewMenu);
	menuBar()->addMenu(helpMenu);
	addToolBar(Qt::TopToolBarArea,toolBar);

	errorPlotWindow = new QWidget(this);
	errorPlotWindow->setWindowFlags(Qt::Window);
	errorPlotWindow->setWindowTitle("Error");

	simPlotWindow = new QWidget(this);
	simPlotWindow->setWindowFlags(Qt::Window);
	simPlotWindow->setWindowTitle("Simulation results");

	fittedPlotWindow = new QWidget(this);
	fittedPlotWindow->setWindowFlags(Qt::Window);
	fittedPlotWindow->setWindowTitle("Data and residuals");

	plotWindowSize = QSize(400,500);

	errorPlot = new DataPlot;
	initialSimPlot = new DataPlot;
	currentSimPlot = new DataPlot;
	residualPlot = new DataPlot;
	dataFitPlot = new DataPlot;

	residualPlot->setTitle("Residuals");
	errorPlot->setTitle("Error");
	dataFitPlot->setTitle("Data and current model");
	initialSimPlot->setTitle("Initial model");
	currentSimPlot->setTitle("Current model");

	QVBoxLayout * plotLayout = new QVBoxLayout;
	plotLayout->addWidget(errorPlot);
	errorPlotWindow->setLayout(plotLayout);
	QCheckBox * checkBox1 = new QCheckBox(" errors ");
	connect(checkBox1,SIGNAL(toggled(bool)),this,SLOT(setErrorPlotVisible(bool)));

	plotLayout = new QVBoxLayout;
	plotLayout->addWidget(initialSimPlot);
	plotLayout->addWidget(currentSimPlot);
	simPlotWindow->setLayout(plotLayout);
	QCheckBox * checkBox2 = new QCheckBox(" model simulation ");
	connect(checkBox2,SIGNAL(toggled(bool)),this,SLOT(setSimPlotVisible(bool)));

	plotLayout = new QVBoxLayout;
	plotLayout->addWidget(dataFitPlot);
	plotLayout->addWidget(residualPlot);
	fittedPlotWindow->setLayout(plotLayout);
	QCheckBox * checkBox3 = new QCheckBox(" data fit ");
	connect(checkBox3,SIGNAL(toggled(bool)),this,SLOT(setFitPlotVisible(bool)));

	QCheckBox * checkBox4 = new QCheckBox(" covariance ");
	covMat = new QTextEdit(this);
	covMat->setReadOnly(true);
	covMat->setWindowFlags(Qt::Window);
	connect(checkBox4,SIGNAL(toggled(bool)),covMat,SLOT(setVisible(bool)));

	paramsTable = new QTableWidget;
	advancedTable = new QTableWidget;
	
	connect(paramsTable,SIGNAL(cellChanged(int,int)),this,SLOT(paramsTableChanged(int,int)));
	connect(advancedTable,SIGNAL(cellChanged(int,int)),this,SLOT(advancedTableChanged(int,int)));
	
	QVBoxLayout * mainLayout = new QVBoxLayout;
	QVBoxLayout * labelsLayout = new QVBoxLayout;
	QVBoxLayout * paramsLayout = new QVBoxLayout;
	
	QGroupBox * labelsGroup = new QGroupBox("");
	labelsLayout->addWidget(sbmlFileLabel = new QLabel("no model loaded"));
	labelsLayout->addWidget(dataFileLabel = new QLabel("no data loaded"));
	labelsGroup->setLayout(labelsLayout);

	QGroupBox * paramsGroup = new QGroupBox(" Parameters ");
	
	QPushButton * randomizeButton = new QPushButton(" Randomize ");
	QDockWidget * randomizerDock = new QDockWidget("Setup randomization");
	
	randomizer = new RandomizerWidget;
	randomizerDock->setWidget(randomizer);
	addDockWidget ( Qt::BottomDockWidgetArea , randomizerDock );
	randomizerDock->setFloating(true);
	randomizerDock->hide();

	QAction * toggleRnd = randomizerDock->toggleViewAction();
	toggleRnd->setChecked(false);
	randomizerDock->hide();
	QToolButton * showHideRnd = new QToolButton;
	showHideRnd->setText(" Configure >>> ");
	showHideRnd->setCheckable(true);
	showHideRnd->setChecked(false);
	viewMenu->addAction(toggleRnd);
	connect(showHideRnd,SIGNAL(toggled(bool)),randomizerDock,SLOT(setVisible(bool)));
	connect(randomizeButton,SIGNAL(clicked()),this,SLOT(randomizeInitialValues()));

	paramsLayout->addWidget(paramsTable);

	QHBoxLayout * randomizeButtonLayout = new QHBoxLayout;
	randomizeButtonLayout->addStretch(3);
	randomizeButtonLayout->addWidget(randomizeButton);
	randomizeButtonLayout->addWidget(showHideRnd);
	randomizeButtonLayout->addStretch(3);
	errorLabel = new QLabel(" error = NA ");
	randomizeButtonLayout->addWidget(errorLabel);
	paramsLayout->addLayout(randomizeButtonLayout);
	paramsGroup->setLayout(paramsLayout);

	QGroupBox * checkboxGroup = new QGroupBox(" Graphs ");
	QGridLayout * checkboxLayout = new QGridLayout;
	checkboxLayout->addWidget(checkBox1,0,0,Qt::AlignLeft);
	checkboxLayout->addWidget(checkBox2,0,1,Qt::AlignLeft);
	checkboxLayout->addWidget(checkBox3,1,0,Qt::AlignLeft);
	checkboxLayout->addWidget(checkBox4,1,1,Qt::AlignLeft);
	checkboxGroup->setLayout(checkboxLayout);

	QGroupBox * algsGroup = new QGroupBox(" Optimization algorithms ");
	QHBoxLayout * algsLayout = new QHBoxLayout;
	algorithmsListBox = new QComboBox;
	algsLayout->addWidget(algorithmsListBox);
	QToolButton * showHideAdv = new QToolButton;
	showHideAdv->setText(" Advanced >>> ");
	showHideAdv->setCheckable(true);
	algsLayout->addWidget(showHideAdv);
	algsLayout->setStretch(0,1);
	algsGroup->setLayout(algsLayout);

	QDockWidget * advancedDock = new QDockWidget("Advanced");
	
	advancedDock->setWidget(advancedTable);
	addDockWidget ( Qt::BottomDockWidgetArea , advancedDock );
	advancedDock->setFloating(true);

	QAction * toggleAdv = advancedDock->toggleViewAction();
	toggleAdv->setChecked(false);
	advancedDock->hide();
	showHideAdv->setChecked(false);
	viewMenu->addAction(toggleAdv);
	connect(showHideAdv,SIGNAL(toggled(bool)),advancedDock,SLOT(setVisible(bool)));

	mainLayout->addWidget(labelsGroup);
	mainLayout->addWidget(paramsGroup);
	mainLayout->addWidget(algsGroup);
	mainLayout->addWidget(checkboxGroup);

	QWidget * widget = new QWidget;
	widget->setLayout(mainLayout);
	setCentralWidget(widget);

	errorPlotWindow->resize(plotWindowSize.width(),plotWindowSize.height()/2);
	simPlotWindow->resize(plotWindowSize);
	fittedPlotWindow->resize(plotWindowSize);
}

void JSimOptimGUI::loadSBML()
{
	QString file(sbmlFile.c_str());
	if (file.isEmpty())
		file = QDir::currentPath();
	else
		file = QString(sbmlFile.c_str());
	
	file = QFileDialog::getOpenFileName(
		this, 
		"select SBML file", 
		file);

	if (file.isNull() || file.isEmpty()) return;

	sbmlFile = std::string(file.toAscii().data());

	connectToRoadRunner();

	if (!roadRunnerLoaded)
		QMessageBox::information(this,"Error","SBW module edu.kgi.roadRunner not loaded");
	else
		sbmlFileLabel->setText(tr("SBML loaded: ") + file);
}

void JSimOptimGUI::loadData()
{
	QString file(dataFile.c_str());
	if (file.isEmpty())
		file = QDir::currentPath();
	else
		file = QString(dataFile.c_str());
	
	file = QFileDialog::getOpenFileName(
		this, 
		"select CSV file", 
		file);

	if (file.isNull() || file.isEmpty()) return;

	QFile reader(file);
	
	if (reader.open(QFile::ReadOnly | QFile::Text))
	{
		varNames.clear();
		DataTable<qreal> data;
		QStringList lines = QString(reader.readAll()).split("\n");
		QStringList words;
		bool ok;
		for (int i=0; i < lines.size(); ++i)
		{
			lines[i].remove('\"');
			words = lines[i].split(",");
			if (words.size() < 1) continue;
			if (i==0)
			{
				data.resize(1,words.size());
				data.setColNames(words);
				for (int j=0; j < words.size(); ++j)
					if (words[j].toLower() != tr("time"))
						varNames.push_back(std::string(words[j].toAscii().data()));
			}
			else
			{
				for (int j=0; j < words.size(); ++j)
				{
					double d = words[j].toDouble(&ok);
					if (ok)
						data.value(i-1,j) = d;
				}
			}
		}
		
		if (data.colName(0).toLower() == tr("time"))
		{
			timeEnd = data.value( data.rows()-1 ,0);
			stepSize = data.value( 1, 0 ) - data.value( 0, 0 );
			timeEnd += stepSize;
		}

		numDataColumns = data.cols()-1;
		dataFitPlot->getDataTable() = data;
		DataTable<qreal> simu_result;

		if (!sbmlFile.empty())
		{
			simu_result = simulate();
			initialSimPlot->plot(simu_result,0,"Initial Simulation");
		}
		createGradedDataPlot(simu_result);
	}

	dataFile = std::string(file.toAscii().data());
	dataFileLabel->setText(tr("Target data: ") + file);
}

void JSimOptimGUI::createGradedDataPlot(DataTable<qreal>& sim)
{
	DataTable<qreal> & data = dataFitPlot->getDataTable();

	QList<QPen> pens;
	int k = 0;
	for (int i=0; i < numDataColumns; ++i,++k)
	{
		if (k >= DataPlot::penList.size())
			k = 0;
		pens << DataPlot::penList[k];
	}

	if (sim.cols() == 0)
	{
		dataFitPlot->plot(data,0,"Target Data",pens);
		return;
	}

	int oldcols = data.cols();
	int newcols = numDataColumns + oldcols;

	DataTable<qreal> residuals;
	residuals.resize(sim.rows(),2);

	QPen pen;
	QColor color;
	data.resize(data.rows(),newcols);

	k = 0;
	while (pens.size() < newcols)
	{
		++k;
		for (int j=0; j < numDataColumns; ++j)
		{
			pen = DataPlot::penList[j];
			color = pen.brush().color();
			color.setAlphaF(0.5*(double)(k*numDataColumns)/(double)newcols);
			pen.setBrush( QBrush(color) );
			pens << pen;
		}
	}

	QStringList colnames(sim.getColNames());
	k = 0;
	for (int i=0; i < numDataColumns; ++i)
	{
		k = colnames.indexOf(data.colName(i+1));
		if (k > -1)
		{
			data.colName(oldcols + i) = colnames[k] + QString::number(iterCount);
			for (int j=0; j < sim.rows() && data.rows(); ++j)
			{
				data.value( j, oldcols + i ) = sim.value(j,k);
				residuals.value( j, 0) = sim.value(j,0);
				residuals.value( j, 1) = data.value(j,i+1) - sim.value(j,k);
			}
		}
	}

	dataFitPlot->plot(data,0,"Target Data",pens);
	pen = QPen(QColor(255,10,10));
	residualPlot->plot(residuals,0,"Residuals",QList<QPen>() << pen);
}

void JSimOptimGUI::run()
{
	++runCount;
	iterCount = 1;

	errorPlot->getDataTable().resize(0,0);
	errorPlot->replot();

	DataTable<qreal> & dat = dataFitPlot->getDataTable();
	dataFitPlot->getDataTable().resize(dat.rows(),numDataColumns+1);

	if (!roadRunnerLoaded)
	{
		QMessageBox::information(this,"Error","RoadRunner module not loaded properly");
	}

	if (!jsimLoaded)
	{
		QMessageBox::information(this,"Error","JSim module not loaded properly");
	}

	if (sbmlFile.empty())
	{
		QMessageBox::information(this,"Error","No SBML file loaded");
	}

	if (dataFile.empty())
	{
		QMessageBox::information(this,"Error","No data file loaded");
	}

	JSimRun * thread = new JSimRun(this);
	thread->jsimModule = jsimOptimizerModule;

	std::vector<std::string> selectedParamNames;
	for (int i=0; i < checkBoxes.size() && i < paramNames.size(); ++i)
		if (checkBoxes[i] && checkBoxes[i]->isChecked())
			selectedParamNames.push_back(paramNames[i]);

	thread->varNames = varNames;
	thread->sbml = sbmlFile;
	thread->csvfile = dataFile;
	thread->paramNames = selectedParamNames;
	thread->varNames = varNames;
	thread->paramInitialValues = initialValues;
	thread->optimParamNames = optimizerParamNames;
	thread->optimParamValues = optimizerParamValues;

	connect(this,SIGNAL(stop()),thread,SLOT(terminate()));

	connect(thread,SIGNAL(errorChanged(double)),this,SLOT(updateErrorPlot(double)));

	connect(thread,SIGNAL(resultsChanged(std::vector<double>)),
			this,SLOT(updateOptimResults(std::vector<double>)));

	thread->start();
}

void JSimOptimGUI::updateErrorPlot(double err)
{
	errorLabel->setText(tr(" iter = ") + QString::number(iterCount) + tr(", error = ") + QString::number(err));
	DataTable<qreal> data = errorPlot->getDataTable();
	data.value(iterCount,0) = iterCount;
	data.value(iterCount,1) = err;
	data.colName(0) = tr("iteration");
	data.colName(runCount) = tr("run ") + QString::number(runCount);
	errorPlot->plot(data,0,"Error");
	++iterCount;

	try
	{
		Service optimService = jsimOptimizerModule.findServiceByName("Optimization");
		ServiceMethod cov = optimService.getMethod("double[][] getCovarianceMatrix()");

		DataBlockReader dataReader = cov.call(DataBlockWriter());
		double** rawData; 
		int numRows=0, numCols=0, i=0;
		dataReader.get(numRows, numCols,rawData);

		QString s;
		for (int i=0; i < numRows; ++i)
		{
			for (int j=0; j < numCols; ++j)
			{
				s += QString::number(rawData[i][j]);
				s += tr("\t");
			}
			s += tr("\n");
		}
		DataBlockReader::free2DArray(numRows, rawData);
		covMat->setText(s);
	}
	catch(SBWException *e)
	{
		std::cout << e->getMessage() << std::endl;
	}
}

void JSimOptimGUI::updateOptimResults(std::vector<double> results)
{
	try
	{
		SBW::connect();
		roadRunnerModule = SBW::getModuleInstance("edu.kgi.roadRunner");		
		Service sim = roadRunnerModule.findServiceByName("sim");
		ServiceMethod loadSBML = sim.getMethod("void loadSBMLFromFile(string)");
		ServiceMethod getParamNames = sim.getMethod("{} getGlobalParameterNames()");
		ServiceMethod getNumParams = sim.getMethod("int getNumberOfGlobalParameters()");
		ServiceMethod setParamValue = sim.getMethod("void setGlobalParameterByIndex(int,double)");
		
		DataBlockWriter emptyArgs;
		QStringList params;

		for (int i=0, j=0; i < checkBoxes.size() && i < paramNames.size(); ++i)
			if (checkBoxes[i] && checkBoxes[i]->isChecked())
			{
				params << QString(paramNames[i].c_str());

				if (initialValues.size() > i && paramsTable->rowCount() > i && results.size() > j)
				{
					initialValues[i] = results[j];
					paramsTable->setItem(i,1, new QTableWidgetItem(QString::number(initialValues[i])));
					++j;
				}
			}

		vector<double> paramValues;
		
		loadSBML.call(DataBlockWriter() << sbmlFile);

		for (int j=0; j < results.size() && j < params.size(); ++j)
			for (int i=0; i < paramNames.size(); ++i)
				if (QString(paramNames[i].c_str()) == params[j])
					setParamValue.call(DataBlockWriter() << i << results[j]);

		DataTable<qreal> simu_result = simulate();
		currentSimPlot->plot(simu_result,0,"Fitted Simulation");
		createGradedDataPlot(simu_result);
	}
	catch(SBWException *e)
	{
		std::cout << e->getMessage() << std::endl;
	}
}

void JSimOptimGUI::paramsTableChanged(int i,int j)
{
	disconnect(paramsTable,SIGNAL(cellChanged(int,int)),this,SLOT(paramsTableChanged(int,int)));

	QTableWidgetItem * cell = 0;
	bool ok = false;
	double d;

	for (int i=0; i < paramsTable->rowCount() && i < paramValues.size(); ++i)
		if (cell = paramsTable->item(i,1))
		{
			QString s = cell->text();
			d = s.toDouble(&ok);
			if (ok)
				paramValues[i] = d;
			else
				cell->setText(QString::number(paramValues[i]));
		}

	connect(paramsTable,SIGNAL(cellChanged(int,int)),this,SLOT(paramsTableChanged(int,int)));
}

void JSimOptimGUI::advancedTableChanged(int i,int j)
{
	disconnect(advancedTable,SIGNAL(cellChanged(int,int)),this,SLOT(advancedTableChanged(int,int)));
	QTableWidgetItem * cell = 0;
	bool ok = false;
	double d;

	optimizerParamNames.resize(advancedTable->rowCount());
	optimizerParamValues.resize(advancedTable->rowCount());

	for (int i=0; i < advancedTable->rowCount(); ++i)
		if (cell = advancedTable->item(i,0))
		{
			QString s = cell->text();
			d = s.toDouble(&ok);
			if (ok)
				optimizerParamValues[i] = d;
			else
				cell->setText(QString::number(optimizerParamValues[i]));
			optimizerParamNames[i] = std::string(advancedTable->verticalHeaderItem(i)->text().toAscii().data());
		}
	connect(advancedTable,SIGNAL(cellChanged(int,int)),this,SLOT(advancedTableChanged(int,int)));
}

void JSimOptimGUI::setErrorPlotVisible(bool b)
{
	errorPlotWindow->setVisible(b);
	if (b)
	{
		//QPoint p = errorPlotWindow->pos();
		errorPlotWindow->resize(plotWindowSize.width(),plotWindowSize.height()/2);
		/*if (simPlotWindow->isVisible())
			p.rx() = simPlotWindow->geometry().right();
		if (fittedPlotWindow->isVisible() &&
			fittedPlotWindow->geometry().right() > p.x())
			p.rx() = fittedPlotWindow->geometry().right();
		errorPlotWindow->move(p - errorPlotWindow->pos());*/
	}
}

void JSimOptimGUI::setSimPlotVisible(bool b)
{
	simPlotWindow->setVisible(b);
	if (b)
	{
		//QPoint p = simPlotWindow->pos();
		simPlotWindow->resize(plotWindowSize);
		/*if (errorPlotWindow->isVisible())
			p.rx() = errorPlotWindow->geometry().right();
		if (fittedPlotWindow->isVisible() &&
			fittedPlotWindow->geometry().right() > p.x())
			p.rx() = fittedPlotWindow->geometry().right();
		simPlotWindow->move(p - simPlotWindow->pos());*/
	}
}

void JSimOptimGUI::setFitPlotVisible(bool b)
{
	fittedPlotWindow->setVisible(b);
	if (b)
	{
		//QPoint p = fittedPlotWindow->pos();
		fittedPlotWindow->resize(plotWindowSize);
		/*if (errorPlotWindow->isVisible())
			p.rx() = errorPlotWindow->geometry().right();
		if (simPlotWindow->isVisible() &&
			simPlotWindow->geometry().right() > p.x())
			p.rx() = fittedPlotWindow->geometry().right();
		fittedPlotWindow->move(p - fittedPlotWindow->pos());*/
	}
}

void JSimOptimGUI::disconnectJsim()
{
	try
	{
		jsimOptimizerModule = SBW::getModuleInstance("JSimSBWOptimServer");
		Service optimService = jsimOptimizerModule.findServiceByName("Optimization");
		ServiceMethod jsimDisconnect = optimService.getMethod("void disconnect()");
		jsimDisconnect.call(DataBlockWriter());
	}
	catch(SBWException *e)
	{
		std::cout << e->getMessage() << std::endl;
	}
}

void JSimOptimGUI::keyPressEvent ( QKeyEvent * event )
{
	emit stop();
	QWidget::keyPressEvent(event);
}

void JSimOptimGUI::mousePressEvent ( QMouseEvent * event )
{
	emit stop();
	QWidget::mousePressEvent(event);
}

