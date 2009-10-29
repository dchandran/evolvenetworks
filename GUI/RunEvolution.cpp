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
	if (args < 5) return 0;
	
	File = argv[1];
	Generations = atoi(argv[2]);
	PopulationSize = atoi(argv[3]);
	StartingPopulationSize = atoi(argv[4]);
	
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
    
	mainWindow.go();
	
	int output = app.exec();
	
    return output;
}

namespace NetworkEvolutionLib
{
	Thread * CurrentThread = 0;
	QSemaphore * semaphore = 0;
	typedef void (*InitFunc)(void);
	
	Thread::Thread(const QString& file, QObject * parent) : QThread(parent)
	{
		lib = new QLibrary(file,this);
		lib->load();
		CurrentThread = this;
	}
	
	void Thread::emitSignal(int iter, int popSz, void** pop, double * fitnessArray, int *** parents)
	{
		if (semaphore)
		{
			semaphore->acquire();
			emit updateScene(iter,popSz,pop,fitnessArray,parents);
			semaphore->acquire();
			semaphore->release();
		}
	}
	
	int MainCallback(int iter, int popSz, GApopulation pop, double * fitnessArray, int *** parents)
	{
		if (CurrentThread)
			CurrentThread->emitSignal(iter,popSz,pop,fitnessArray,parents);
		if (callback)
			return callback(iter,popSz,pop,fitnessArray,parents);
		return 0;
	}
	
	void Thread::run()
	{
		if (!CurrentThread || !lib || !lib->isLoaded()) return;
		
		void * f0 = lib->resolve("init");
		void * f1 = lib->resolve("fitness");
		void * f2 = lib->resolve("callback");
		if (f0 && f1 && f2)
		{
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
	
	QSize MainWindow::sizeHint() const
	{
		return QSize(400,400);
	}

	MainWindow::MainWindow()
	{
		if (semaphore)
			delete semaphore;
		
		semaphore = new QSemaphore(1);
		QGraphicsView * view = new GraphicsView(scene = new QGraphicsScene, this);
		scene->setBackgroundBrush(QBrush(Qt::black));
		view->setDragMode(QGraphicsView::ScrollHandDrag);
		setCentralWidget(view);
	}
	
	MainWindow::~MainWindow()
	{
		if (semaphore)
			delete semaphore;
		QList<QGraphicsItem*> items = scene->items();
		for (int i=0; i < items.size(); ++i)
			if (items[i] && !items[i]->parentItem())
				delete items[i];
	}
	
	void MainWindow::insertTextItem(double number, int line, int rank, int parent1, int parent2)
	{
		double w = 10.0, h = 10.0;
		//QGraphicsSimpleTextItem * textItem = new QGraphicsSimpleTextItem(QString::number(number));
		
		//textItem->setBrush(QBrush(QColor(200,255,200)));
		//textItem->setPos( 10.0 + rank * 100.0, 10.0 + line * 100.0 );
		
		QGraphicsRectItem * rectItem = new QGraphicsRectItem(	100.0 + rank * w * 2.0, 100.0 + line * h * 5.0,
																w, h );
																//2.0 * textItem->boundingRect().width(), 
																//2.0 * textItem->boundingRect().height());
		
		//textItem->scale(2.0,2.0);
		
		rectItem->setBrush(QBrush(QColor((int)(255 * number),10,(int)((1.0-number)*255))));
		rectItem->setPen(QPen(QColor(100,255,100)));
		
		scene->addItem(rectItem);
		//scene->addItem(textItem);
		
		if (parent1 > 0)
		{
			QGraphicsLineItem * lineItem = new QGraphicsLineItem(
									100.0 + parent1 * w * 2.0, 100.0 + h + (line-1) * h * 5.0,
									100.0 + rank * w * 2.0, 100.0 - h + (line) * h * 5.0);
									
			lineItem->setPen(QPen(QColor(100,255,100)));
			scene->addItem(lineItem);
		}
		
		if (parent2 > 0)
		{
			QGraphicsLineItem * lineItem = new QGraphicsLineItem(
									100.0 + parent2 * w * 2.0, 100.0 + h + (line-1) * h * 5.0,
									100.0 + rank * w * 2.0, 100.0 - h + (line) * h * 5.0);
			lineItem->setPen(QPen(QColor(255,100,100)));
			scene->addItem(lineItem);
		}
	}
	
	void MainWindow::updateScene(int iter, int popSz, GApopulation P, double * scores, int *** parents)
	{
		if (!fitness) return;
		
		double min, max;
		
		for (int i=0; i < popSz; ++i)
		{
			if (i==0)
				min = max = scores[i];
			else
			{
				if (scores[i] > max) max = scores[i];
				if (scores[i] < min) min = scores[i];
			}
		}
		
		max -= min;
		
		for (int i=0; i < popSz; ++i)
		{
			int * p = getImmediateParents(i,iter, parents);
			if (p && p[0])
			{
				insertTextItem((scores[i]-min)/max, iter, i, p[0], p[1]);
			}
			else
			{
				insertTextItem((scores[i]-min)/max, iter, i, 0, 0);
			}
			if (p)
				delete p;
		}
		
		delete scores;
		
		if (semaphore)
			semaphore->release();
	}
	
	void MainWindow::go()
	{
		show();
		
		if (!File) return;
		
		Thread * thread = new Thread( tr(File) , this);
		connect(thread,SIGNAL(updateScene(int , int, void** , double *, int *** )),this,SLOT(updateScene(int , int, void** , double *, int ***)));
		thread->start();
	}	
}
