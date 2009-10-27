#include "RunEvolution.h"
#include "reactionNetwork.h"

using namespace NetworkEvolutionLib;
using namespace Tinkercell;

char * File = 0;
int Generations = 0;
int PopulaionSize = 0;
int StartingPopulaionSize = 0;
GACallbackFnc callback = 0;

int main(int args, char *argv[])
{
	if (args < 4) return 0;
	
	File = argv[0];
	Generations = atoi(argv[1]);
	PopulaionSize = atoi(argv[2]);
	StartingPopulaionSize = atoi(argv[3]);
	
    QApplication app(argc, argv);
	
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
	
	void MainWindow::MainCallback
	
	QSize MainWindow::sizeHint() const
	{
		return QSize(400,400);
	}

	MainWindow::MainWindow()
	{
		QGraphicsView * view = new QGraphicsView(scene = new QGraphicsScene, this);
		scene->setBackgroundBrush(QBrush(Qt::black));
		view->setDragMode(QGraphicsView::ScrollHandDrag);
		setCentralWidget(view);
		
		insertTextItem(5.0,1,5,0,0);
		
		insertTextItem(12.0,2,1,5,0);
		insertTextItem(10.0,2,2,5,0);
		insertTextItem(1.0,2,3,5,0);
		insertTextItem(42.0,2,4,5,0);
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
	
	typedef (*InitFunc)(void);
	
	void MainWindow::go()
	{
		if (!File) return;
		
		QString filename(File);
		void * f0 = QLibrary::resolve(filename, "init");
		void * f1 = QLibrary::resolve(filename, "fitness");
		void * f2 = QLibrary::resolve(filename, "callback");
		if (f0 && f1 && f2)
		{
			InitFunc initFunc = (InitFunc)f0;
			GAFitnessFnc fitnessFunc = (GAFitnessFnc)f1;
			callback = (GACallbackFnc)f2;
				
			GApopulation P;
			initFunc();
			P = evolveNetworks(StartingPopulationSize,PopulationSize,Generations,&fitnessFunc,&MainCallback);
			GAfree(P);
		}
	}	
}
