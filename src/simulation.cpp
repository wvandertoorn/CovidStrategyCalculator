#include "simulation.h"
#include "mainwindow.h"

#include <unsupported/Eigen/MatrixFunctions>

#include <QVBoxLayout>
#include <QHeaderView>

#include <QtCharts/QAreaSeries>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QLegendMarker>
#include <QtCharts/QValueAxis>

std::vector<int> Simulation::sub_compartments = {5,1,13,1,1}; //predetect, presympt, symp, post, risk-node (sink is dropped)
int Simulation::nr_compartments = std::accumulate(sub_compartments.begin(), sub_compartments.end(), 0);
float Simulation::t_inf = 100.;

Simulation::Simulation(MainWindow *parent) : QObject(parent)
{
    this->m_parent = parent;
    collect_data(this->m_parent);

    initial_states.setZero(nr_compartments);

    switch (this->mode) {
        case 0: initial_states(0) = 1.0; // mode exposure
                break;
        case 1:
        {
                int first_symptomatic_compartment = sub_compartments[0] + sub_compartments[1];
                this->initial_states(first_symptomatic_compartment) = 1.0; // mode symptom onset
                initial_states(nr_compartments - 1) = 1.0;
                //TODO read only location?
                break;
        }
        case 2: break;
    }
}

Simulation::Simulation(MainWindow *parent, Eigen::MatrixXf init_states, float prerisk) : QObject(parent)
{
    this->m_parent = parent;
    collect_data(this->m_parent);

    // Eigen::MatrixXf tmp;
    // tmp.setZero(nr_compartments, 1);
    init_states.resize(nr_compartments, 1);
    // tmp(Eigen::seq(0, Eigen::last-1), 0).array() = init_states.array();
    init_states(nr_compartments - 1, 0) = init_states(Eigen::seq(sub_compartments[0],
                                                                 sub_compartments[0]
                                                                 + sub_compartments[1]
                                                                 + sub_compartments[2] -1), 0)
                                                                .sum();
    //TODO read only location?

    this->initial_states = init_states;

    this->pre_test_infect_prob = prerisk;
}

std::tuple<std::vector<Eigen::MatrixXf>,
           std::vector<Eigen::MatrixXf>> Simulation::run_prevalence(std::vector<float> week_incidences)
{
    time_passed = 0;
    this->fraction_asymtomatic = 1.;
    //TODO document the fact that the prevalence estimator does not assume symptomatic screening in the first simulation.

    Eigen::MatrixXf states_today_mean, states_today_lev, states_today_uev,
                    probs_today_mean, probs_today_lev, probs_today_uev;

    states_today_mean.setZero(1, nr_compartments);
    states_today_lev.setZero(1, nr_compartments);
    states_today_uev.setZero(1, nr_compartments);

    float total_residence_mean = std::accumulate(residence_times_mean.begin(), residence_times_mean.end(), 0);
    // float tot_lev = std::accumulate(residence_times_lev.begin(), residence_times_lev.end(), 0);
    // float tot_uev = std::accumulate(residence_times_uev.begin(), residence_times_uev.end(), 0);

    for (int week = 0; week < 5; ++week)
    {
        float day_incidence = week_incidences[week] / 7;
        this->quarantine = (week + 1) * 7 - 1;

        // set intial states according to incidence report
        int counter = 0;
        for (int i = 0; i < 4; ++i)
        {
            for (int j=0; j < sub_compartments[i]; ++j)
            {
                    float probability = residence_times_mean[i] / total_residence_mean * day_incidence / sub_compartments[i];
                    initial_states(counter) = probability;
                    counter++;
            }
        }

        this->run();
        states_today_mean.array() += X_mean(Eigen::seq(Eigen::last-6, Eigen::last), Eigen::all).colwise().sum().array();
        states_today_lev.array() += X_lev(Eigen::seq(Eigen::last-6, Eigen::last), Eigen::all).colwise().sum().array();
        states_today_uev.array() += X_uev(Eigen::seq(Eigen::last-6, Eigen::last), Eigen::all).colwise().sum().array();
        //TODO think about how to handle risk node in prevalence estimator

        initial_states.array() = 0.0; //needed to set initial state of sink compartment back to zero
    }

    probs_today_mean = assemble_phases(states_today_mean, sub_compartments);
    probs_today_lev = assemble_phases(states_today_lev, sub_compartments);
    probs_today_uev = assemble_phases(states_today_uev, sub_compartments);

    std::vector<Eigen::MatrixXf> states{states_today_mean, states_today_lev, states_today_uev};
    std::vector<Eigen::MatrixXf> probs{probs_today_mean, probs_today_lev, probs_today_uev};

    return std::make_tuple(states, probs);
}

void Simulation::collect_data(MainWindow *parent)
{
    // strategy
    mode = parent->mode_ComboBox->currentIndex();
    mode_str = parent->mode_ComboBox->currentText().toStdString();
    time_passed = parent->time_passed->value();
    quarantine = parent->quarantine->value();
    use_symptomatic_screening = parent->use_symptomatic_screening->isChecked();

    // parameters
    // mean prediction
    residence_times_mean.clear();
    residence_times_mean.push_back(parent->percentage_predetection->value()/100 * parent->inc_mean->value());
    residence_times_mean.push_back((1- parent->percentage_predetection->value()/100) * parent->inc_mean->value());
    residence_times_mean.push_back(parent->symp_mean->value());
    residence_times_mean.push_back(parent->post_mean->value());

    // lower extreme value prediction
    residence_times_lev.clear();
    residence_times_lev.push_back(parent->percentage_predetection->value()/100 * parent->inc_lev->value());
    residence_times_lev.push_back((1- parent->percentage_predetection->value()/100) * parent->inc_lev->value());
    residence_times_lev.push_back(parent->symp_lev->value());
    residence_times_lev.push_back(parent->post_mean->value());

    // upper extreme value prediction
    residence_times_uev.clear();
    residence_times_uev.push_back(parent->percentage_predetection->value()/100 * parent->inc_uev->value());
    residence_times_uev.push_back((1- parent->percentage_predetection->value()/100) * parent->inc_uev->value());
    residence_times_uev.push_back(parent->symp_uev->value());
    residence_times_uev.push_back(parent->post_mean->value());

    fraction_asymtomatic = parent->percentage_asymptomatic->value() /100;

    sensitivity = parent->pcr_sens->value() /100;
    specificity = parent->pcr_spec->value() /100;

    this->test_type = parent->test_type->currentText();

    switch (parent->test_type->currentIndex()) {
        case 1: sensitivity = sensitivity * (parent->rel_antigen_sens->value() /100);
                break; //antigen
    }

    if (use_symptomatic_screening){
        fraction_asymtomatic = parent->percentage_asymptomatic->value() /100;
    } else {fraction_asymtomatic = 1.;}

    t_test = collect_t_test(parent->test_date_checkboxes);
}

std::vector<int> Simulation::collect_t_test(std::vector<QCheckBox*> boxes)
{
    std::vector<int> v{};
    int counter = time_passed;
    for (auto box : boxes)
    {
        if ( box->isChecked() )
            v.push_back(counter);
        ++counter;
    }
    return v;
}

Eigen::VectorXf Simulation::calc_rates(std::vector<float> times,
                                       std::vector<int> comp
                                       )
{
    int n = std::accumulate(comp.begin(), comp.end(), 0);
    Eigen::VectorXf rates(n-1);
    int counter = 0;
    for (int i=0; i < 4; ++i)
    {
        for (int j=0; j < comp[i]; ++j)
        {
                rates[counter] = comp[i] /times[i];
                counter++;
        }
    }
    return rates;
}

Eigen::MatrixXf Simulation::calc_S(int n)
{
    Eigen::MatrixXf S;
    S.setZero(n, n-1);

    for (int i=0; i < n-1; ++i)
    {
        S(i, i) = -1;
        S(i+1, i) = 1;
    }
    return S;
}

Eigen::MatrixXf Simulation::calc_A(Eigen::MatrixXf S,
                                   Eigen::VectorXf r,
                                   std::vector<int> comp)
{
    int n = S.rows();
    Eigen::MatrixXf temp;
    temp.setZero(n, n-1);

    for (int i = 0; i < n;  ++i)
    {
        for (int j=0; j < n-1; ++j)
        {
            temp(i,j) = r(j);
        }
    }

    Eigen::MatrixXf A;
    A = S.cwiseProduct(temp);

    Eigen::MatrixXf A_square(n-1, n-1);
    A_square = A(Eigen::seq(0, Eigen::last-1), Eigen::all); //drop last row, sink node not of interest

    Eigen::MatrixXf A_augmented;
    A_augmented.setZero(n,n);
    A_augmented(Eigen::seq(0, Eigen::last-1),Eigen::seq(0, Eigen::last-1)) = A_square;
    A_augmented(Eigen::last, Eigen::seq(comp[0], comp[0] + comp[1] - 1)).array() = 1.;
    A_augmented(Eigen::last, Eigen::seq(comp[0] + comp[1], comp[0] + comp[1] + comp[2] - 1)).array() = fraction_asymtomatic;
    return A_augmented;
}

Eigen::MatrixXf  Simulation::calc_X(int delay,
                                    int qrntn,
                                    Eigen::MatrixXf A,
                                    Eigen::VectorXf states
                                    )
{
    int t_end = delay + qrntn + 1;

    Eigen::VectorXf time(t_end);
    for (int i=0; i< t_end; ++i)
    {
        time(i) = float(i);
    }

    int n = states.size();

    Eigen::MatrixXf X;
    X.setZero(t_end, n);

    for (int i=0; i<t_end; ++i)
    {
        X.row(i) =(A * time[i]).exp()*states;
    }

    return X;
}

float Simulation::calc_risk_at_T(Eigen::MatrixXf A,
                                Eigen::VectorXf states,
                                float time_T)
{
    int n = states.size();
    Eigen::MatrixXf X;
    X.setZero(1, n);

    X.row(0) = (A * time_T).exp()*states;

    return X(0, Eigen::last);
}

Eigen::MatrixXf Simulation::assemble_phases(Eigen::MatrixXf X,
                                            std::vector<int> comp
                                            )
{
    int n_time = X.rows();
    Eigen::MatrixXf assembled(n_time, 5);
    int col_counter=0;
    for (int i=0; i<5; ++i)
    {
        assembled(Eigen::all, i) = X(Eigen::all, Eigen::seq(col_counter,
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

    QtCharts::QLineSeries *risk_low = new QtCharts::QLineSeries;
    QtCharts::QLineSeries *risk_mean = new QtCharts::QLineSeries;
    QtCharts::QLineSeries *risk_high = new QtCharts::QLineSeries;
    risk_mean->setName("Is- or will become infectious");

    for (int j=0; j<n_time; ++j)
    {
        float m = result_matrix_mean(j, Eigen::seq(0,2)).sum();
        float lev = result_matrix_lev(j, Eigen::seq(0,2)).sum();
        float uev = result_matrix_uev(j, Eigen::seq(0,2)).sum();

        risk_low->append(j-time_passed, lev);
        risk_mean->append(j-time_passed, m);
        risk_high->append(j-time_passed, uev);

    }
    QtCharts::QAreaSeries *risk_area = new QtCharts::QAreaSeries(risk_low, risk_high);

    QColor risk_color(255, 128, 128, 128);
    QPen risk_pen(risk_color);
    risk_area->setPen(risk_pen);
    QBrush risk_brush(risk_color);
    risk_area->setBrush(risk_brush);
    QPen red_pen(Qt::red);
    risk_mean->setPen(red_pen);

    chart->addSeries(risk_area);
    chart->addSeries(risk_mean);
    chart->legend()->markers(risk_area)[0]->setVisible(false);

    //------------------

    QtCharts::QLineSeries *detect_low = new QtCharts::QLineSeries;
    QtCharts::QLineSeries *detect_high = new QtCharts::QLineSeries;

    Eigen::MatrixXf detectable = calculate_assay_sensitivity();
    for (int j=0; j<n_time; ++j)
    {
        detect_low->append(j-time_passed, detectable(j,0));
        detect_high->append(j-time_passed, detectable(j,1));
    }
    QtCharts::QAreaSeries *detect_area = new QtCharts::QAreaSeries(detect_low, detect_high);
    detect_area->setName("Infected and detectable");

    QColor detect_color(140, 129, 152, 128);
    QPen detect_pen(detect_color);
    detect_area->setPen(detect_pen);
    QBrush detect_brush(detect_color);
    detect_area->setBrush(detect_brush);

    chart->addSeries(detect_area);

    chart->createDefaultAxes();
    chart->axes(Qt::Horizontal).first()->setTitleText("Day");
    chart->axes(Qt::Vertical).first()->setTitleText("Probability");

    QtCharts::QValueAxis *axisX = qobject_cast<QtCharts::QValueAxis*>(chart->axes(Qt::Horizontal).first());
    axisX->setTickCount(n_time);

    chart->legend()->show();
    QFont f("Helvetica", 10);
    chart->legend()->setFont(f);

    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    return chartView;
}

QTableWidget* Simulation::create_table()
{
    int n_time = result_matrix_mean.rows();
    int offset = time_passed;

    QStringList h_labels;
    for (int i = - offset; i < n_time - offset; ++i)
    {
        h_labels << QString::number(i);
    }

    QTableWidget *table = new QTableWidget(1, n_time);
    table->setHorizontalHeaderLabels(h_labels);
    table->setVerticalHeaderLabels((QStringList() << "Test efficacy"));

    float start_mean = result_matrix_mean(0, Eigen::seq(0, Eigen::last-1)).sum();
    float start_lev = result_matrix_lev(0, Eigen::seq(0, Eigen::last-1)).sum();
    float start_uev = result_matrix_uev(0, Eigen::seq(0, Eigen::last-1)).sum();

    for (int i = 0; i< n_time; ++i)
    {
        float perc_mean = ((1-specificity) * result_matrix_mean(i, 0) +
                          (result_matrix_mean(i, Eigen::seq(1, 3)).sum() * sensitivity)) / start_mean * 100;
        float perc_lev = ((1-specificity) * result_matrix_lev(i, 0) +
                          (result_matrix_lev(i,  Eigen::seq(1, 3)).sum() * sensitivity)) / start_lev * 100;
        float perc_uev = ((1-specificity) * result_matrix_uev(i, 0) +
                          (result_matrix_uev(i, Eigen::seq(1, 3)).sum() * sensitivity)) / start_uev * 100;
        table->setItem(0, i, new QTableWidgetItem(QString::number(perc_mean, 'f', 2) + "%"
                                                  + "  (" + QString::number(perc_uev, 'f', 2)
                                                  + ", " + QString::number(perc_lev, 'f', 2) + ")"));

        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    return table;
}

void Simulation::create_result_log()
{
    QLabel *label = new QLabel(tr("Result log"));

    QTableWidget *table = new QTableWidget(1, 9);
    table->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    table->setHorizontalHeaderLabels((QStringList() << "mode"
                                                    << "time passed \n [days]"
                                                    << "quarantine/isolation \n [days]"
                                                    << "sympt screening \n (perc asympt)"
                                                    << "test \n [days]"
                                                    << "test type"
                                                    << "pre-procedure \n risk [%]"
                                                    << "residual risk \n [%]"
                                                    << "fold risk \n reduction"));

    write_row_result_log(table);
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    connect(
        table->horizontalHeader(),
        SIGNAL(sectionResized(int, int, int)),
        table,
        SLOT(resizeRowsToContents()));

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(label, 0);
    layout->addWidget(table, 90);

    this->m_parent->log->setLayout(layout);
}

void Simulation::write_row_result_log(QTableWidget *table)
{
    table->setItem(0, 0, new QTableWidgetItem(QString(this->mode_str.c_str())));
    table->setItem(0, 1, new QTableWidgetItem(QString::number(this->time_passed)));
    table->setItem(0, 2, new QTableWidgetItem(QString::number(this->quarantine)));

    QString boolText = this->use_symptomatic_screening ? "yes" : "no";

    if (this->use_symptomatic_screening){
        table->setItem(0, 3, new QTableWidgetItem(boolText
                                                 + " ("
                                                 + QString::number(this->fraction_asymtomatic * 100)
                                                 + "%)"));
    } else {table->setItem(0, 3, new QTableWidgetItem(boolText));}

    if (this->t_test.size())
    {
        QString days{};
        for (int day : this->t_test)
        {
            days += (QString::number(day - this->time_passed) + ", ");
        }
        days.chop(2);
        table->setItem(0, 4, new QTableWidgetItem(days));

        table->setItem(0, 5, new QTableWidgetItem(this->test_type));
    }
    else
    {
      table->setItem(0, 4, new QTableWidgetItem());
      table->setItem(0, 5, new QTableWidgetItem());
    }

    table->setItem(0, 6, new QTableWidgetItem(QString::number(pre_test_infect_prob *100, 'f', 2)));

    table->setItem(0, 7, new QTableWidgetItem(QString::number(calculate_strategy_result(result_matrix_mean)* 100,
                                                             'f', 2)
                                             + "  ("
                                             + QString::number(calculate_strategy_result(result_matrix_uev)* 100,
                                                               'f', 2)
                                             + ", "
                                             + QString::number(calculate_strategy_result(result_matrix_lev)* 100,
                                                               'f', 2)
                                             + ")" ));
    table->setItem(0, 8, new QTableWidgetItem(QString::number(pre_test_infect_prob
                                                               / calculate_strategy_result(result_matrix_mean),
                                                             'f', 2)
                                             + "  ("
                                             + QString::number(pre_test_infect_prob
                                                                 / calculate_strategy_result(result_matrix_uev),
                                                               'f', 2)
                                             + ", "
                                             + QString::number(pre_test_infect_prob
                                                                 / calculate_strategy_result(result_matrix_lev),
                                                               'f', 2)
                                             + ")" ));

    for (int i = 0; i < 9; ++i)
    {
        table->item(0,i)->setFlags(table->item(0,i)->flags() &  ~Qt::ItemIsEditable);
    }
}

void Simulation::update_result_log()
{
    QWidget *widget = this->m_parent->log->layout()->takeAt(1)->widget();
    QTableWidget* table = qobject_cast<QTableWidget*>(widget);
    table->insertRow(0);
    write_row_result_log(table);

    this->m_parent->log->layout()->addWidget(table);
}

float Simulation::calculate_strategy_result(Eigen::MatrixXf matrix)
{
    float risk = 1.;

    if (t_test.size())
    {
        for (int day : t_test)
        {
            risk = risk * (1. - sensitivity) * (matrix(day, 1) + fraction_asymtomatic * matrix(day, 2) ) +
                   risk * specificity * matrix(day, 0);
        }
        if (t_test.back() < matrix.rows()-1)
        {
            risk = risk * (matrix(Eigen::last, Eigen::seq(0, 1)).sum() + fraction_asymtomatic * matrix(Eigen::last, 2))
                    / (matrix(t_test.back(), Eigen::seq(0, 1)).sum() + fraction_asymtomatic * matrix(t_test.back(), 2));
        }
    }
    else risk = matrix(Eigen::last, Eigen::seq(0, 1)).sum() + fraction_asymtomatic * matrix(Eigen::last, 2);
    return risk;
}

Eigen::MatrixXf Simulation::calculate_assay_sensitivity()
{
    int n_time = result_matrix_mean.rows();
    Eigen::MatrixXf detectibility(n_time, 2);

    for (int j=0; j<n_time; ++j)
    {
        float m = result_matrix_mean(j, 0) * (1 - specificity) + result_matrix_mean(j, Eigen::seq(1, 3)).sum() * sensitivity;
        float lev = assay_detectibility_worst_case(j, 0) * (1 - specificity) + assay_detectibility_worst_case(j, Eigen::seq(1, 3)).sum() * sensitivity;
        float uev = assay_detectibility_best_case(j, 0) * (1 - specificity) + assay_detectibility_best_case(j, Eigen::seq(1, 3)).sum() * sensitivity;

        std::vector<float> v{m, lev, uev};

        int maxElementIndex = std::max_element(v.begin(), v.end()) - v.begin();
        int minElementIndex = std::min_element(v.begin(), v.end()) - v.begin();

        detectibility(j, 0) =  v[minElementIndex];
        detectibility(j, 1) =  v[maxElementIndex];
    }
    return detectibility;
}

void Simulation::output_results()
{
    QtCharts::QChartView* plot = create_plot();
    QTableWidget* table = create_table();

    QVBoxLayout *layout = new QVBoxLayout();
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setSpacing( 0 );

    layout->addWidget(plot, 85);
    layout->addWidget(table, 8);

    // if previous layout exists, delete
    if (this->m_parent->chart->layout())
    {
        QLayout *hb = this->m_parent->chart->layout();
        while(!hb->isEmpty())
        {
            QWidget *w = hb->takeAt(0)->widget();
            delete w;
        }
        delete hb;
    }

    this->m_parent->chart->setLayout(layout);
    this->m_parent->chart->update();

    // initialize or update
    if (!this->m_parent->log->layout())
    {
        create_result_log();
    }
    else
    {
        update_result_log();
    }

}

void Simulation::run()
{
    Eigen::VectorXf rates;
    Eigen::MatrixXf S;
    Eigen::MatrixXf A;
    Eigen::MatrixXf X;

    rates = calc_rates(residence_times_mean, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates, sub_compartments);
    this->X_mean = calc_X(time_passed, quarantine, A, initial_states);
    this->result_matrix_mean = assemble_phases(X_mean, sub_compartments);
    this->t_inf_risk_mean = calc_risk_at_T(A, initial_states, t_inf);
    this->pre_procedure_risk_mean = t_inf_risk_mean - initial_states(Eigen::last, 0);

    rates = calc_rates(residence_times_lev, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates, sub_compartments);
    this->X_lev = calc_X(time_passed, quarantine, A, initial_states);
    this->result_matrix_lev = assemble_phases(X_lev, sub_compartments);
    this->t_inf_risk_lev = calc_risk_at_T(A, initial_states, t_inf);
    this->pre_procedure_risk_lev = t_inf_risk_lev - initial_states(Eigen::last, 0);

    rates = calc_rates(residence_times_uev, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates, sub_compartments);
    this->X_uev = calc_X(time_passed, quarantine, A, initial_states);
    this->result_matrix_uev = assemble_phases(X_uev, sub_compartments);
    this->t_inf_risk_uev = calc_risk_at_T(A, initial_states, t_inf);
    this->pre_procedure_risk_uev = t_inf_risk_uev - initial_states(Eigen::last, 0);

    /* time dependent sensitivity assay requires (up, low, low, low) and (low, up, up, up)
     *  for worst and best case respectively.
     */
    std::vector<float> worst_case_extreme{residence_times_uev[0],
                                          residence_times_lev[1],
                                          residence_times_lev[2],
                                          residence_times_lev[3]};
    std::vector<float> best_case_extreme{residence_times_lev[0],
                                         residence_times_uev[1],
                                         residence_times_uev[2],
                                         residence_times_uev[3]};

    rates = calc_rates(worst_case_extreme, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates, sub_compartments);
    X = calc_X(time_passed, quarantine, A, initial_states);
    this->assay_detectibility_worst_case = assemble_phases(X, sub_compartments);

    rates = calc_rates(best_case_extreme, sub_compartments);
    S = calc_S(nr_compartments);
    A = calc_A(S, rates, sub_compartments);
    X = calc_X(time_passed, quarantine, A, initial_states);
    this->assay_detectibility_best_case = assemble_phases(X, sub_compartments);

    calculate_assay_sensitivity();
}
