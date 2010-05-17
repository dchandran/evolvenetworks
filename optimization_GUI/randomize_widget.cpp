#include <QLabel>
#include "randomize_widget.h"

double RandomizerWidget::rnorm()
{
	double r = ((double)rand())/(double)RAND_MAX;
	return (sqrt(-2.0*log( r )) *cos (2.0*3.14159* r ));
}

RandomizerWidget::RandomizerWidget()
{
	QVBoxLayout * layout = new QVBoxLayout;

	comboBox = new QComboBox;
	connect(comboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(currentIndexChanged(int)));

	QHBoxLayout * comboBoxLayout = new QHBoxLayout;
	QGroupBox * group1 = new QGroupBox(" Select parameter ");
	comboBoxLayout->addWidget(comboBox);
	group1->setLayout(comboBoxLayout);

	useNorm = new QRadioButton("Normal random number");
	useUnif = new QRadioButton("Uniform random number");
	multiplyCurrent = new QCheckBox(" Use multiples of current value? ");

	connect(useNorm,SIGNAL(toggled(bool)),this,SLOT(radioToggled(bool)));

	meanLine = new QLineEdit("");
	connect(meanLine,SIGNAL(textChanged(const QString&)),this,SLOT(muChanged(const QString&)));

	sdLine = new QLineEdit("");
	connect(sdLine,SIGNAL(textChanged(const QString&)),this,SLOT(sdChanged(const QString&)));

	minLine = new QLineEdit("");
	connect(minLine,SIGNAL(textChanged(const QString&)),this,SLOT(minChanged(const QString&)));

	maxLine = new QLineEdit("");
	connect(maxLine,SIGNAL(textChanged(const QString&)),this,SLOT(maxChanged(const QString&)));

	QGridLayout * normunifLayout = new QGridLayout;
	QGroupBox * group2 = new QGroupBox(" Random variable options ");
	normunifLayout->addWidget(useNorm,0,0);
	normunifLayout->addWidget(useUnif,0,1);
	normunifLayout->addWidget(multiplyCurrent,1,0);
	group2->setLayout(normunifLayout);

	QGridLayout * linesLayout = new QGridLayout;
	QGroupBox * group3 = new QGroupBox(" Random variable options ");
	linesLayout->addWidget(new QLabel("mean:"),0,0);
	linesLayout->addWidget(meanLine,0,1);
	linesLayout->addWidget(new QLabel("st.dev:"),0,2);
	linesLayout->addWidget(sdLine,0,3);
	linesLayout->addWidget(new QLabel("min:"),1,0);
	linesLayout->addWidget(minLine,1,1);
	linesLayout->addWidget(new QLabel("max:"),1,2);
	linesLayout->addWidget(maxLine,1,3);
	group3->setLayout(linesLayout);

	layout->addWidget(group1);
	layout->addWidget(group2);
	layout->addWidget(group3);
	setLayout(layout);
}

void RandomizerWidget::initialize(const QStringList & list, const QList<double>& values)
{
	//typeOfRandom.clear();
	names.clear();
	means.clear();
	sds.clear();
	mins.clear();
	maxs.clear();

	if (list.size() != values.size()) return;

	names = list;
	means = values;
	for (int i=0; i < names.size(); ++i)
	{
		mins << values[i]/10.0;
		maxs << values[i]*10.0;
		sds << values[i]*2.0;
		//typeOfRandom << 1;
	}

	multiplyCurrent->setChecked(true);
	comboBox->addItems(list);

	if (list.size() > 0)
		comboBox->setCurrentIndex(0);
}

double RandomizerWidget::getRandomValue(const QString& name, double currentValue)
{
	int i = names.indexOf(name);
	return getRandomValue(i,currentValue);
}

double RandomizerWidget::getRandomValue(int i, double currentValue)
{	
	if (i < 0) return currentValue;

	if (useNorm->isChecked())
	{
		if (multiplyCurrent->isChecked())
			return (currentValue * rnorm()*sds[i]);
		else
			return (means[i] + rnorm()*sds[i]);
	}
	else
	{
		double r = ((double)rand())/(double)RAND_MAX;
		if (multiplyCurrent->isChecked())
			return (currentValue * r * (maxs[i]-mins[i]));
		else
			return (mins[i] + r * (maxs[i]-mins[i]));
	}

	return currentValue;
}

void RandomizerWidget::radioToggled(bool)
{
	//if (typeOfRandom.size() < 1) return;

	int i = comboBox->currentIndex();
	/*if (i >= typeOfRandom.size())
		return;*/

	if (!useNorm->isChecked())
	{
		//typeOfRandom[i] = 0;
		meanLine->setEnabled(false);
		sdLine->setEnabled(false);
		minLine->setEnabled(true);
		maxLine->setEnabled(true);
	}
	else
	{
		//typeOfRandom[i] = 1;
		meanLine->setEnabled(true);
		sdLine->setEnabled(true);
		minLine->setEnabled(false);
		maxLine->setEnabled(false);
	}
}

void RandomizerWidget::currentIndexChanged(int i)
{
	/*if (typeOfRandom.size() > i)
	{
		useNorm->setChecked(typeOfRandom[i]==0);
		useUnif->setChecked(typeOfRandom[i]==1);
	}*/

	if (means.size() > i)
		meanLine->setText(QString::number(means[i]));

	if (sds.size() > i)
		sdLine->setText(QString::number(sds[i]));

	if (mins.size() > i)
		minLine->setText(QString::number(mins[i]));

	if (maxs.size() > i)
		maxLine->setText(QString::number(maxs[i]));
}

void RandomizerWidget::muChanged(const QString& s)
{
	bool ok;
	double d = s.toDouble(&ok);
	if (!ok) return;

	int i = comboBox->currentIndex();
	if (i >= means.size())
		return;

	means[i] = d;
}

void RandomizerWidget::sdChanged(const QString& s)
{
	bool ok;
	double d = s.toDouble(&ok);
	if (!ok) return;

	int i = comboBox->currentIndex();
	if (i >= sds.size())
		return;

	sds[i] = d;
}

void RandomizerWidget::minChanged(const QString& s)
{
	bool ok;
	double d = s.toDouble(&ok);
	if (!ok) return;

	int i = comboBox->currentIndex();
	if (i >= mins.size())
		return;

	mins[i] = d;
}

void RandomizerWidget::maxChanged(const QString& s)
{
	bool ok;
	double d = s.toDouble(&ok);
	if (!ok) return;

	int i = comboBox->currentIndex();
	if (i >= maxs.size())
		return;

	maxs[i] = d;
}
