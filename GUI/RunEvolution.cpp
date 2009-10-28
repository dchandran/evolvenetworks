#include <iostream>
#include "RunEvolution.h"
extern "C"
{
#include "reactionNetwork.h"
}

using namespace std;
using namespace NetworkEvolutionLib;
using namespace Tinkercell;

char * File = 0;
int Generations = 0;
int PopulationSize = 0;
int StartingPopulationSize = 0;
GACallbackFnc callback = 0;
GAFitnessFnc fitness = 0;

int main(int args, char *argv[])
{
	if (args < 4) return 0;
	
	File = argv[0];
	Generations = atoi(argv[1]);
	PopulationSize = atoi(argv[2]);
	StartingPopulationSize = atoi(argv[3]);
	
    QApplication app(args, argv);
	
	QString appDir = QCoreApplication::applicationDirPath();
	
    QFile styleFile(appDir + QString("/networkevolve.qss"));
	
	if (styleFile.open(QFile::ReadOnly | QFile::Text))
    {
        app.setStyleSheet(styleFile.readAll());
        styleFile.close();
    }

    MainWindow mainWindow;

    mainWindow.show();
    
	int output = app.exec();
	
    return output;
}

namespace NetworkEvolutionLib
{
	MainWindow * MainWindow::mainWindow = 0;

	int MainWindow::MainCallback(int iter, GApopulation pop, int popSz)
	{
		if (mainWindow)
			mainWindow->updateScene(iter,pop,popSz);
	}
	
	QSize MainWindow::sizeHint() const
	{
		return QSize(400,400);
	}

	MainWindow::MainWindow()
	{
		mainWindow = this;
		
		QGraphicsView * view = new QGraphicsView(scene = new QGraphicsScene, this);
		scene->setBackgroundBrush(QBrush(Qt::black));
		view->setDragMode(QGraphicsView::ScrollHandDrag);
		setCentralWidget(view);
		
		go();
	}
	
	MainWindow::~MainWindow()
	{
		QList<QGraphicsItem*> items = scene->items();
		for (int i=0; i < items.size(); ++i)
			if (items[i] && !items[i]->parentItem())
				delete items[i];
	}
	
	void MainWindow::insertTextItem(double number, int line, int rank, int parent1, int parent2)
	{
		QGraphicsSimpleTextItem * textItem = new QGraphicsSimpleTextItem(QString::number(number));
		
		textItem->setBrush(QBrush(QColor(200,255,200)));
		textItem->setPos( 10.0 + rank * 100.0, 10.0 + line * 100.0 );
		
		QGraphicsRectItem * rectItem = new QGraphicsRectItem(	10.0 + rank * 100.0, 10.0 + line * 100.0,
																2.0 * textItem->boundingRect().width(), 
																2.0 * textItem->boundingRect().height());
		
		textItem->scale(2.0,2.0);
		
		rectItem->setBrush(QBrush(QColor(100,100,100)));
		rectItem->setPen(QPen(QColor(100,255,100)));
		
		scene->addItem(rectItem);
		scene->addItem(textItem);
		
		if (parent1 > 0)
		{
			QGraphicsLineItem * lineItem = new QGraphicsLineItem(
									10.0 + parent1 * 100.0, 30.0 + (line-1) * 100.0,
									10.0 + rank * 100.0, 10.0 + (line) * 100.0);
			lineItem->setPen(QPen(QColor(100,255,100)));
			scene->addItem(lineItem);
		}
		
		if (parent2 > 0)
		{
			QGraphicsLineItem * lineItem = new QGraphicsLineItem(
									10.0 + parent2 * 100.0, 30.0 + (line-1) * 100.0,
									10.0 + rank * 100.0, 10.0 + (line) * 100.0);
			lineItem->setPen(QPen(QColor(255,100,100)));
			scene->addItem(lineItem);
		}
	}
	
	void MainWindow::updateScene(int iter, GApopulation P, int popSz)
	{
		if (!fitness) return;
		
		for (int i=0; i < popSz; ++i)
		{
			double score = fitness(P[i]);
			int * p = getImmediateParents(iter,i);
			if (p && p[0])
			{
				insertTextItem(score, iter, i, p[0], p[1]);
			}
			else
			{
				insertTextItem(score, iter, i, 0, 0);
			}
		}
	}
	
	typedef void (*InitFunc)(void);
	
	void MainWindow::go()
	{
		if (!File) return;
		
		cout << "go started\n";
		
		QString filename(File);
		
		QLibrary * lib = new QLibrary(tr(File),this);
		
		if (!lib->isLoaded())
		{
			delete lib;
			return;
		}
		
		void * f0 = lib->resolve(filename, "init");
		void * f1 = lib->resolve(filename, "fitness");
		void * f2 = lib->resolve(filename, "callback");
		if (f0 && f1 && f2)
		{
			cout << "all three f's\n";
			InitFunc initFunc = (InitFunc)f0;
			fitness = (GAFitnessFnc)f1;
			callback = (GACallbackFnc)f2;
				
			GApopulation P;
			initFunc();
			P = evolveNetworks(StartingPopulationSize,PopulationSize,Generations,fitness,&MainCallback);
			GAfree(P);
		}
		
		fitness = 0;
		callback = 0;
		delete lib;
	}	
}
