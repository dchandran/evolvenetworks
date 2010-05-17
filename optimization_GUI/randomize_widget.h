#ifndef RANDOMIZE_WIDGET
#define RANDOMIZE_WIDGET

#include <math.h>
#include <QWidget>
#include <QPushButton>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QList>
#include <QStringList>
#include <QComboBox>
#include <QLineEdit>
#include <QGridLayout>
#include <QRadioButton>
#include <QCheckBox>

class RandomizerWidget : public QWidget
{
	Q_OBJECT

public:

	double rnorm();
	RandomizerWidget();
	void initialize(const QStringList & list, const QList<double>& values);
	double getRandomValue(const QString& name, double currentValue);
	double getRandomValue(int i, double currentValue);

public slots:

	void radioToggled(bool);
	void currentIndexChanged(int i);
	void muChanged(const QString& s);
	void sdChanged(const QString& s);
	void minChanged(const QString& s);
	void maxChanged(const QString& s);

private:
	QRadioButton * useNorm;
	QRadioButton * useUnif;
	QCheckBox * multiplyCurrent;
	QComboBox * comboBox;
	QLineEdit * meanLine;
	QLineEdit * sdLine;
	QLineEdit * minLine;
	QLineEdit * maxLine;

	QStringList names;
	//QList<int> typeOfRandom;
	QList<double> means, sds, mins, maxs;
};
#endif
