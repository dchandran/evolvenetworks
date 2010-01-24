#include <iostream>
#include "RunEvolution.h"
extern "C"
{
#include "blocks.h"
}

using namespace std;
using namespace NetworkEvolutionLib;
using namespace Tinkercell;

char * File = 0;
int Generations = 0;
int PopulationSize = 0;
int StartingPopulationSize = 0;
double MaxFitness = 0.0;
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
	typedef double (*MaxFitnessFunction)(void);
	
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
		
		MaxFitness = 0.0;
		fitness = 0;
		callback = 0;
		
		void * f0 = lib->resolve("init");
		void * f1 = lib->resolve("fitness");
		void * f2 = lib->resolve("callback");
		void * f3 = lib->resolve("maxfitness");
		if (f0 && f1)
		{
			InitFunc initFunc = (InitFunc)f0;
			fitness = (GAFitnessFnc)f1;
			
			if (f3)
			{
				MaxFitnessFunction f = (MaxFitnessFunction)(f3);
				MaxFitness = f();
			}
			
			if (f2)
			{
				callback = (GACallbackFnc)f2;
			}
			
			GApopulation P;
			initFunc();
			P = evolveNetworks(StartingPopulationSize,PopulationSize,Generations,fitness,&MainCallback);
			GAfree(P);
		}
		
		MaxFitness = 0.0;
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
	
	void MainWindow::insertTextItem(double fitness, double max, int iteration, int index, int parent1, int parent2)
	{
		double w = 10.0, h = 10.0;
		
		if (index < 0 || index >= nextGen.size()) return;
		
		QColor color((int)(255 * fitness/max), (int)((1.0-fitness/max)*255), 0);
		
		if ( iteration == 0 )
		{
			double x = 1000.0 * mtrand(), y = 1000.0 * mtrand();
			QGraphicsRectItem * rectItem = new QGraphicsRectItem(x, y, w, h);
			rectItem->setBrush(QBrush(color));
			rectItem->setPen(Qt::NoPen);
			scene->addItem(rectItem);
			rectItem->setToolTip(QString::number(fitness));
			
			nextGen[index] = QPointF(x,y);
			return;
		}
		
		if (parent1 >= 0 && parent1 < previousGen.size())
		{
			double x = previousGen[parent1].x() + 3.0 * w * cos(2 * 3.14159 * mtrand()), 
				   y = previousGen[parent1].y() + 3.0 * h * sin(2 * 3.14159 * mtrand());
			QGraphicsRectItem * rectItem = new QGraphicsRectItem(x,y, w, h);
			rectItem->setBrush(QBrush(color));
			rectItem->setPen(Qt::NoPen);
			scene->addItem(rectItem);
			rectItem->setToolTip(QString::number(fitness));
			
			nextGen[index] = QPointF(x,y);
		}
		
		if (parent2 >= 0 && parent2 < previousGen.size())
		{
			double x = previousGen[parent2].x() + 3.0 * w * cos(2 * 3.14159 * mtrand()), 
				   y = previousGen[parent2].y() + 3.0 * h * sin(2 * 3.14159 * mtrand());
			QGraphicsRectItem * rectItem = new QGraphicsRectItem(x,y, w, h);
			rectItem->setBrush(QBrush(color));
			rectItem->setPen(Qt::NoPen);
			scene->addItem(rectItem);
			rectItem->setToolTip(QString::number(fitness));
			
			nextGen[index] = QPointF(x,y);
		}
	}
	
	void MainWindow::updateScene(int iter, int popSz, GApopulation P, double * scores, int *** parents)
	{
		double max = scores[0];
		
		if (MaxFitness > 0.0) max = MaxFitness;
		
		nextGen.clear();
		
		for (int i=0; i < popSz; ++i)
			nextGen << QPointF(1000.0 * mtrand(), 1000.0 * mtrand());
		
		int sz = popSz/2;
		
		if (iter > 0)
			sz = popSz/4;
		
		for (int i=0; i < sz; ++i)
		{
			int * p = getImmediateParents(i,iter, parents);
			if (p && p[0])
			{
				insertTextItem(scores[i], max, iter, i, p[0]-1, p[1]-1);
			}
			else
			{
				insertTextItem(scores[i], max, iter, i, -1, -1);
			}
			if (p)
			{
				delete p;
			}
		}
		
		previousGen = nextGen;
		
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
