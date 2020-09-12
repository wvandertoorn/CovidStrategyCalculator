#include "simulation.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <QVBoxLayout>
#include <QTableWidget>
#include <QDesktopWidget>
#include <QHeaderView>

#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

std::vector<int> Simulation::sub_compartments = {1,1,11,1,1};
int Simulation::nr_compartments = 15;

Simulation::Simulation() : QObject()
{
    initial_states.setZero(nr_compartments);
    initial_states(0) = 1.0;
}

Simulation::Simulation(MainWindow *parent) : QObject(parent)
{
    initial_states.setZero(nr_compartments);
    initial_states(0) = 1.0;

    m_parent = parent;
    collect_data(m_parent);
}

void Simulation::collect_data(MainWindow *parent)
{
    // input
    time_passed_known = parent->ui->checkBox_known->isChecked();
    if (time_passed_known){
        time_passed = parent->ui->input_timepassed->value();
    } else {time_passed = 0;}
    pre_test_infect_prob = 1;
    quarantaine = parent->ui->input_quarantaine->value();

    // parameters
    residence_times.clear();
    residence_times.push_back(parent->ui->param_time_predetect->value());
    residence_times.push_back(parent->ui->param_time_presymp->value());
    residence_times.push_back(parent->ui->param_time_infect->value());
    residence_times.push_back(parent->ui->param_time_post->value());

    percentage_asymt = parent->ui->param_perc_asym->value() /100;
    pcr_sensitivity = parent->ui->param_pcr_sens->value() /100;
    pcr_specificity = parent->ui->param_pcr_spec->value() /100;

    t_test = collect_t_test(parent->test_boxes);
}

std::vector<int> Simulation::collect_t_test(std::vector<QCheckBox*> boxes){
    std::vector<int> v{};
    int counter = time_passed;
    for (auto box : boxes){
        if ( box->isChecked() )
            v.push_back(counter);
        ++counter;
    }
    return v;
}

Eigen::VectorXf Simulation::calc_rates(std::vector<float> times,
                                       std::vector<int> comp
                                       ){
    int n = std::accumulate(comp.begin(), comp.end(), 0);
    Eigen::VectorXf r(n-1);
    int counter = 0;
    for (int i=0; i < 4; ++i){
        for (int j=0; j < comp[i]; ++j){
                r[counter] = comp[i] /times[i];
                counter++;
        }
    }
    return r;
}

Eigen::MatrixXf Simulation::calc_S(int n){

    Eigen::MatrixXf m;
    m.setZero(n, n-1);

    for (int i=0; i < n-1; ++i){
        m(i, i) = -1;
        m(i+1, i) = 1;
    }
    return m;
}

Eigen::MatrixXf Simulation::calc_A(Eigen::MatrixXf S_,
                                   Eigen::VectorXf r){
    int n = S_.rows();
    Eigen::MatrixXf temp;
    temp.setZero(n, n-1);

    for (int i = 0; i < n;  ++i){
        for (int j=0; j < n-1; ++j){
            temp(i,j) = r(j);
        }
    }

    return S_.cwiseProduct(temp);
}

Eigen::MatrixXf  Simulation::calc_X(float delay,
                                    float qrntn,
                                    Eigen::MatrixXf A_,
                                    Eigen::VectorXf states
                                    ){
    int t_end = std::ceil(delay) + std::floor(qrntn) + 1;

    Eigen::VectorXf time(t_end);
    for (int i=0; i< t_end; ++i){
        time(i) = float(i);
    }

    int n = states.size();

    Eigen::MatrixXf A_square(n-1,n-1);
    A_square = A_(Eigen::seq(0, Eigen::last-1), Eigen::all); //drop last row

    Eigen::MatrixXf x_analytical;
    x_analytical.setZero(t_end, n-1);

    for (int i=0; i<t_end; ++i){
        x_analytical.row(i) =(A_square * time[i]).exp()*states.head(n-1);
    }

    return x_analytical;
}

Eigen::MatrixXf Simulation::assemble_phases(Eigen::MatrixXf X_,
                                            std::vector<int> comp
                                            ){
    int n_time = X_.rows();
    Eigen::MatrixXf assembled(n_time, 4);
    int col_counter=0;
    for (int i=0; i<4; ++i){

        assembled(Eigen::all, i) = X_(Eigen::all, Eigen::seq(col_counter,
                                                             col_counter + comp.at(i) -1))
                                    .rowwise()
                                    .sum();
        col_counter = col_counter + comp.at(i);
    }
    return assembled;
}

QtCharts::QChartView* Simulation::create_plot()
{
    int n_time = result_matrix.rows();

    std::vector<QString> species_names{"incubation",
                                       "pre-(a)symptomatic",
                                       "symptomatic",
                                       "post-symptomatic"};
    std::vector<QColor> species_colors{QColor(173, 245, 66, 255),
                                       QColor(245, 188, 66, 255),
                                       QColor(66, 179, 245, 255),
                                       QColor(255, 77, 169, 255)};

    QtCharts::QLineSeries *asymptomatic = new QtCharts::QLineSeries();
    asymptomatic->setName(tr("asymptomatic"));
    asymptomatic->setPen(Qt::DashLine);
    asymptomatic->setColor(QColor(30, 0, 227, 255));

    QtCharts::QChart *chart = new QtCharts::QChart();

    for (int i=0; i<4; ++i){
        QtCharts::QLineSeries *series = new QtCharts::QLineSeries();
        series->setName(species_names.at(i));
        series->setPen(Qt::DashLine);
        series->setColor(species_colors.at(i));

        if (i==2){
            // symptomatic, split/scale by proportion asymptomatic
            for (int j=0; j<n_time; ++j){
                series->append(j - time_passed, (1- percentage_asymt)* result_matrix(j, i));
                asymptomatic->append(j - time_passed, (percentage_asymt)* result_matrix(j, i));
            }
            chart->addSeries(series);
            chart->addSeries(asymptomatic);
        }
        else {
            for (int j=0; j<n_time; ++j){
                series->append(j - time_passed, result_matrix(j, i));
            }
            chart->addSeries(series);
        }
    }

    // plot categories
    QPen fat_pen;
    fat_pen.setWidth(3);
    QtCharts::QLineSeries *potential = new QtCharts::QLineSeries();
    potential->setName("Is or will become infectious");

    potential->setPen(fat_pen);
    potential->setColor("black");
    for (int j=0; j<n_time; ++j){
        potential->append(j - time_passed, result_matrix(j, Eigen::seq(0,2)).sum());
    }
    chart->addSeries(potential);

    QtCharts::QLineSeries *infectious = new QtCharts::QLineSeries();
    infectious->setName("Infectious & detectable");
    infectious->setPen(fat_pen);
    infectious->setColor("red");
    for (int j=0; j<n_time; ++j){
        infectious->append(j - time_passed, result_matrix(j, Eigen::seq(1,2)).sum()*pcr_sensitivity);
    }
    chart->addSeries(infectious);

    chart->createDefaultAxes();
    QtCharts::QValueAxis *axisX = new QtCharts::QValueAxis;
    axisX->setRange(-std::ceil(time_passed), quarantaine);
    int tick_count = std::ceil(time_passed) + quarantaine + 1;
    axisX->setTickCount(tick_count);
    axisX->setMinorTickCount(1);
    axisX->setLabelFormat(tr("%i"));
    axisX->setTitleText(tr("Day"));
    chart->setAxisX(axisX);

    for (auto axis : chart->axes(Qt::Vertical)){
        axis->setTitleText("Probability");
    }

    chart->legend()->show();
    QFont f("Helvetica", 10);
    chart->legend()->setFont(f);

    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    return chartView;
}

QTableWidget* Simulation::create_table(){
    int n_time = result_matrix.rows();
    int offset = time_passed;

    QStringList h_labels;
    for (int i = - offset; i< n_time - offset; ++i){
        h_labels << QString::number(i);
    }

    QTableWidget *table = new QTableWidget(1, n_time);
    table->setHorizontalHeaderLabels(h_labels);
    table->setVerticalHeaderLabels((QStringList() << "Infectious and detectable"));

    for (int i = 0; i< n_time; ++i){
        table->setItem(0, i, new QTableWidgetItem(QString::number(result_matrix(i, Eigen::seq(1,2)).sum()
                                                                  * pcr_sensitivity * 100, 'f', 2) + "%"));


        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
    if (n_time < 20){
        table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    }
    return table;
}

void Simulation::create_result_log(){
    QLabel *label = new QLabel(tr("Result log"));

    QTableWidget *table = new QTableWidget(1, 5);
    table->setHorizontalHeaderLabels((QStringList() << "time passed (days)"
                                                    << "quarantaine (days)"
                                                    << "test moments (days)"
                                                    << "residual risk (%)"
                                                    << "risk reduction (factor)"));

    write_row_result_log(table);
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(label,0);
    layout->addWidget(table,90);

    m_parent->ui->result_log->setLayout(layout);
}
void Simulation::write_row_result_log(QTableWidget *table){
    table->setItem(0,0, new QTableWidgetItem(QString::number(time_passed)));
    table->setItem(0,1, new QTableWidgetItem(QString::number(quarantaine)));
    if (t_test.size()){
        QString days{};
        for (int day : t_test){
            days += (QString::number(day-time_passed) + ", ");
        }
        days.chop(2);
        table->setItem(0,2, new QTableWidgetItem(days));
    }
    else {
        table->setItem(0,2, new QTableWidgetItem());
    }

    float result = calculate_strategy_result();
    table->setItem(0,3, new QTableWidgetItem(QString::number(result*100, 'f', 2)));
    table->setItem(0,4, new QTableWidgetItem(QString::number(pre_test_infect_prob/result, 'f', 2)));

    for (int i=0; i<5; ++i){
        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
}
void Simulation::update_result_log(){
    QWidget *widget = m_parent->ui->result_log->layout()->takeAt(1)->widget();
    QTableWidget* table = qobject_cast<QTableWidget*>(widget);
    table->insertRow(0);
    write_row_result_log(table);

    m_parent->ui->result_log->layout()->addWidget(table);
}

float Simulation::calculate_strategy_result(){
    float risk = 1.;

    if (t_test.size()){
        for (int day : t_test){
            risk = risk * (1.-pcr_sensitivity) * result_matrix(day, Eigen::seq(1,2)).sum() +
                   risk * pcr_specificity * result_matrix(day, 0);
        }
        if (t_test.back() < result_matrix.rows()-1){
            risk = risk * result_matrix(Eigen::last, Eigen::seq(0,2)).sum() /result_matrix(t_test.back(), Eigen::seq(0,2)).sum();
        }
    }
    else risk = result_matrix(Eigen::last, Eigen::seq(0,2)).sum();
    return risk;
}

void Simulation::output_results(){
    QtCharts::QChartView* plot = create_plot();
    QTableWidget* table = create_table();

    QVBoxLayout *layout = new QVBoxLayout();
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setSpacing( 0 );

    layout->addWidget(plot, 85);
    layout->addWidget(table, 8);

    // if previous layout exists, delete
    if (m_parent->ui->results->layout()){
        QLayout *hb = m_parent->ui->results->layout();
        while(!hb->isEmpty()) {
            QWidget *w = hb->takeAt(0)->widget();
            delete w;
        }
        delete hb;
    }
    else {
        m_parent->ui->residual_risk->setText("Residual risk:");
        m_parent->ui->risk_reduction->setText("Risk reduction:");
    }
    m_parent->ui->results->setLayout(layout);
    m_parent->ui->results->update();

    float result = calculate_strategy_result();
    m_parent->ui->result_residual_risk->setText(QString::number(result*100, 'f', 2) + "%");
    m_parent->ui->result_risk_reduction->setText(QString::number(pre_test_infect_prob/result, 'f', 2));

    if (!m_parent->ui->result_log->layout()){
        create_result_log();
    }
    else {
        update_result_log();
    }

}

void Simulation::run(){
    rates = calc_rates(residence_times, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates);
    X = calc_X(time_passed, quarantaine, A, initial_states);
    result_matrix = assemble_phases(X, sub_compartments);
    output_results();
}
