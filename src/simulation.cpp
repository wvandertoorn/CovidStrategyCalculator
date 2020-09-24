#include "simulation.h"
#include "mainwindow.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <QVBoxLayout>
#include <QHeaderView>

#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QCandlestickSeries>
#include <QtCharts/QCandlestickSet>
#include <QtCharts/QLineSeries>
#include <QtCharts/QLegendMarker>
#include <QtCharts/QValueAxis>
#include <QtCharts/QBarCategoryAxis>

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
    time_passed = parent->time_passed->value();
    pre_test_infect_prob = 1;
    quarantine = parent->quarantine->value();

    // parameters
    residence_times_mean.clear();
    residence_times_mean.push_back(parent->pred_mean->value()/100 * parent->inc_mean->value());
    residence_times_mean.push_back((1- parent->pred_mean->value()/100) * parent->inc_mean->value());
    residence_times_mean.push_back(parent->symp_mean->value());
    residence_times_mean.push_back(parent->post_mean->value());

    residence_times_lev.clear();
    residence_times_lev.push_back(parent->pred_mean->value()/100 * parent->inc_lev->value());
    residence_times_lev.push_back((1- parent->pred_mean->value()/100) * parent->inc_lev->value());
    residence_times_lev.push_back(parent->symp_lev->value());
    residence_times_lev.push_back(parent->post_lev->value());

    residence_times_uev.clear();
    residence_times_uev.push_back(parent->pred_mean->value()/100 * parent->inc_uev->value());
    residence_times_uev.push_back((1- parent->pred_mean->value()/100) * parent->inc_uev->value());
    residence_times_uev.push_back(parent->symp_uev->value());
    residence_times_uev.push_back(parent->post_uev->value());

    percentage_asymt = parent->asymp_mean->value() /100;
    pcr_sensitivity = parent->pcr_sens->value() /100;
    pcr_specificity = parent->pcr_spec->value() /100;

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
    int n_time = result_matrix_mean.rows();
    QtCharts::QChart *chart = new QtCharts::QChart();

    QtCharts::QCandlestickSeries *risk_error = new QtCharts::QCandlestickSeries();
    risk_error->setName("Is or will become infectious");
    QPen red_pen(Qt::red,3);
    QBrush red_brush(Qt::red);
    risk_error->setPen(red_pen);
    risk_error->setBrush(red_brush);
    for (int j=0; j<n_time; ++j){
        float lev = result_matrix_lev(j, Eigen::seq(0,2)).sum();
        float uev = result_matrix_uev(j, Eigen::seq(0,2)).sum();
        float m = result_matrix_mean(j, Eigen::seq(0,2)).sum();

        QtCharts::QCandlestickSet *set = new QtCharts::QCandlestickSet(m,uev,lev,m, j-time_passed);
        risk_error->append(set);
    }
    chart->addSeries(risk_error);

    QtCharts::QCandlestickSeries *detectable_error = new QtCharts::QCandlestickSeries();
    detectable_error->setName("Detectable");
    QPen black_pen(Qt::black,3);
    QBrush black_brush(Qt::black);
    detectable_error->setPen(black_pen);
    detectable_error->setBrush(black_brush);
    for (int j=0; j<n_time; ++j){
        float lev = result_matrix_lev(j, Eigen::seq(1,2)).sum()*pcr_sensitivity;
        float uev = result_matrix_uev(j, Eigen::seq(1,2)).sum()*pcr_sensitivity;
        float m = result_matrix_mean(j, Eigen::seq(1,2)).sum()*pcr_sensitivity;

        QtCharts::QCandlestickSet *set = new QtCharts::QCandlestickSet(m,uev,lev,m, j-time_passed);
        detectable_error->append(set);
    }
    chart->addSeries(detectable_error);

    QStringList categories;
    for (int j=0; j<n_time; ++j){
        categories << QString::number(j-time_passed);
    }
    QtCharts::QBarCategoryAxis *axisX_cat = new QtCharts::QBarCategoryAxis;
    axisX_cat->setCategories(categories);
    chart->addAxis(axisX_cat, Qt::AlignBottom);
    axisX_cat->setTitleText("Day");
    risk_error->attachAxis(axisX_cat);
    detectable_error->attachAxis(axisX_cat);

    QtCharts::QValueAxis *axisY = new QtCharts::QValueAxis();
    chart->addAxis(axisY, Qt::AlignLeft);
    risk_error->attachAxis(axisY);
    detectable_error->attachAxis(axisY);
    axisY->setRange(0, 1);
    axisY->setTitleText("Probability");

    chart->legend()->show();
    QFont f("Helvetica", 10);
    chart->legend()->setFont(f);

    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    return chartView;
}

QTableWidget* Simulation::create_table(){
    int n_time = result_matrix_mean.rows();
    int offset = time_passed;

    QStringList h_labels;
    for (int i = - offset; i< n_time - offset; ++i){
        h_labels << QString::number(i);
    }

    QTableWidget *table = new QTableWidget(1, n_time);
    table->setHorizontalHeaderLabels(h_labels);
    table->setVerticalHeaderLabels((QStringList() << "Risk posing and detectable"));

    for (int i = 0; i< n_time; ++i){
        float perc_mean = (1-pcr_specificity) * result_matrix_mean(i, 0) +
                          result_matrix_mean(i, Eigen::seq(1,2)).sum()
                          * pcr_sensitivity * 100;
        float perc_lev = (1-pcr_specificity) * result_matrix_lev(i, 0) +
                          result_matrix_lev(i, Eigen::seq(1,2)).sum()
                          * pcr_sensitivity * 100;
        float perc_uev = (1-pcr_specificity) * result_matrix_uev(i, 0) +
                          result_matrix_uev(i, Eigen::seq(1,2)).sum()
                          * pcr_sensitivity * 100;
        table->setItem(0, i, new QTableWidgetItem(QString::number(perc_mean, 'f', 2) + "%"
                                                  + "  (" + QString::number(perc_uev, 'f', 2)
                                                  + ", " + QString::number(perc_lev, 'f', 2) + ")"));


        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    return table;
}

void Simulation::create_result_log(){
    QLabel *label = new QLabel(tr("Result log"));

    QTableWidget *table = new QTableWidget(1, 5);
    table->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    table->setHorizontalHeaderLabels((QStringList() << "time passed"
                                                    << "quarantine"
                                                    << "test"
                                                    << "risk reduction [%]"
                                                    << "risk reduction [factor]"));

    write_row_result_log(table);
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    connect(
        table->horizontalHeader(),
        SIGNAL(sectionResized(int, int, int)),
        table,
        SLOT(resizeRowsToContents()));

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(label,0);
    layout->addWidget(table,90);

    m_parent->log->setLayout(layout);
}

void Simulation::write_row_result_log(QTableWidget *table){
    table->setItem(0,0, new QTableWidgetItem(QString::number(time_passed)));
    table->setItem(0,1, new QTableWidgetItem(QString::number(quarantine)));
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

    table->setItem(0,3, new QTableWidgetItem(QString::number((pre_test_infect_prob-calculate_strategy_result(result_matrix_mean))/pre_test_infect_prob*100, 'f', 2)
                                             + "  (" + QString::number((pre_test_infect_prob-calculate_strategy_result(result_matrix_uev))/pre_test_infect_prob*100, 'f', 2)
                                             + ", " + QString::number((pre_test_infect_prob-calculate_strategy_result(result_matrix_lev))/pre_test_infect_prob*100, 'f', 2)
                                             + ")" ));
    table->setItem(0,4, new QTableWidgetItem(QString::number(pre_test_infect_prob/calculate_strategy_result(result_matrix_mean), 'f', 2)
                                             + "  (" + QString::number(pre_test_infect_prob/calculate_strategy_result(result_matrix_uev), 'f', 2)
                                             + ", " + QString::number(pre_test_infect_prob/calculate_strategy_result(result_matrix_lev), 'f', 2)
                                             + ")" ));

    for (int i=0; i<5; ++i){
        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
}

void Simulation::update_result_log(){
    QWidget *widget = m_parent->log->layout()->takeAt(1)->widget();
    QTableWidget* table = qobject_cast<QTableWidget*>(widget);
    table->insertRow(0);
    write_row_result_log(table);

    m_parent->log->layout()->addWidget(table);
}

float Simulation::calculate_strategy_result(Eigen::MatrixXf matrix){
    float risk = 1.;

    if (t_test.size()){
        for (int day : t_test){
            risk = risk * (1.-pcr_sensitivity) * matrix(day, Eigen::seq(1,2)).sum() +
                   risk * pcr_specificity * matrix(day, 0);
        }
        if (t_test.back() < matrix.rows()-1){
            risk = risk * matrix(Eigen::last, Eigen::seq(0,2)).sum() /matrix(t_test.back(), Eigen::seq(0,2)).sum();
        }
    }
    else risk = matrix(Eigen::last, Eigen::seq(0,2)).sum();
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
    if (m_parent->chart->layout()){
        QLayout *hb = m_parent->chart->layout();
        while(!hb->isEmpty()) {
            QWidget *w = hb->takeAt(0)->widget();
            delete w;
        }
        delete hb;
    }
    m_parent->chart->setLayout(layout);
    m_parent->chart->update();

    if (!m_parent->log->layout()){
        create_result_log();
    }
    else {
        update_result_log();
    }

}

void Simulation::run(){
    Eigen::VectorXf rates;
    Eigen::MatrixXf S;
    Eigen::MatrixXf A;
    Eigen::MatrixXf X;

    rates = calc_rates(residence_times_mean, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates);
    X = calc_X(time_passed, quarantine, A, initial_states);
    result_matrix_mean = assemble_phases(X, sub_compartments);

    rates = calc_rates(residence_times_lev, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates);
    X = calc_X(time_passed, quarantine, A, initial_states);
    result_matrix_lev = assemble_phases(X, sub_compartments);

    rates = calc_rates(residence_times_uev, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates);
    X = calc_X(time_passed, quarantine, A, initial_states);
    result_matrix_uev = assemble_phases(X, sub_compartments);

    output_results();
}
