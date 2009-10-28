#ifndef NETWORK_EVOLVE_VISUALIZE_GUI_H
#define NETWORK_EVOLVE_VISUALIZE_GUI_H

#include <QApplication>
#include <QCoreApplication>
#include <QWidget>
#include <QMainWindow>
#include <QSettings>
#include <QLabel>
#include <QTextEdit>
#include <QPushButton>
#include <QSplitter>
#include <QString>
#include <QStringList>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QGroupBox>
#include <QComboBox>
#include <QLineEdit>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QHeaderView>
#include <QToolBar>
#include <QFile>
#include <QProcess>
#include <QTextStream>
#include <QFileInfo>
#include <QDir>
#include <QFileDialog>
#include <QLibrary>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QGraphicsSimpleTextItem>
#include <QGraphicsRectItem>
#include <QGraphicsLineItem>
#include "CodeEditor.h"

namespace NetworkEvolutionLib
{
	class MainWindow : public QMainWindow
	{
		Q_OBJECT
		
	public:
		MainWindow();
		~MainWindow();
		QSize sizeHint() const;
	
		static int MainCallback(int, void**, int);

	public slots:
		void go();
		void updateScene(int , void** , int );
		
	private:
		
		void insertTextItem(double, int, int, int, int);
		QGraphicsScene * scene;
		QList<QGraphicsItem*> previousLine;
		
	};
}

#endif

