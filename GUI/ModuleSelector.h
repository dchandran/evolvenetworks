/****************************************************************************
** 
**
****************************************************************************/

#ifndef NETWORK_EVOLUTION_GUI_MODULE_SELECTOR_H
#define NETWORK_EVOLUTION_GUI_MODULE_SELECTOR_H

#include <QWidget>
#include <QTableWidget>
#include <QTabWidget>
#include <QToolButton>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLyaout>


#ifdef Q_WS_WIN
#define MY_EXPORT __declspec(dllexport)
#else
#define MY_EXPORT
#endif

namespace NetworkEvolutionLib
{
	class MY_EXPORT ModuleSelector : public QWidget
	{
		Q_OBJECT
	public:
		ModuleSelector(QWidget * parent = 0);
		
	};

}

#endif

