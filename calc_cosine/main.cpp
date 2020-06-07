#include <QCoreApplication>
#include <QtCore>


//Class to store basic information about an analysis
class Analysis {
public:
    Analysis();
    Analysis(QString name_) {name = name_;}
    void set_metadata(QString key, QString value);

    QMap<int,double> buckets;
    QMap<QString, QString> metadata;


    QString name;
    float ecol_act;
    float msme_act;
    float saur_act;
    QString strain;
    QString media;
    QString supp;
    QString sfi_id;

};

//Helper function to read data from Bruker ProfileAnalysis output
QList<Analysis> read_file_profile_analysis(QString name) {
    QList<Analysis> ret;
    QFile file(name);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return ret;
    bool first=true;
    while(!file.atEnd()) {
        QByteArray line = file.readLine();
        if(first) {
            first=false;
            continue;
        }
        QStringList splitted = QString(line).split(";");
        //qDebug() << "Num: " << splitted.size();
        QString name = splitted.at(0);
        name = name.replace("\"","");
        Analysis ana(name);
        for(int i=1; i<splitted.size(); i++) {
            ana.buckets[i-1] = splitted.at(i).toDouble();
        }
        ret.append(ana);
    }
    return ret;
}

//Helper function to read data from Bruker MetaboScape output
QList<Analysis> read_file_metaboscape(QString name) {
    QList<Analysis> ret;
    QFile file(name);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return ret;

    while(!file.atEnd()) {

        QByteArray line = file.readLine();
        QStringList splitted = QString(line).split(",");
        //qDebug() << "Num: " << splitted.size();
        QString name = splitted.at(0);
        name = name.replace("\"","");
        Analysis ana(name);
        for(int i=1; i<splitted.size(); i++) {
            ana.buckets[i-1] = splitted.at(i).toDouble();
        }
        ret.append(ana);
    }
    return ret;
}

//Helper function to read data from XCMS output
QList<Analysis> read_file_xcms(QString name) {
    QList<Analysis> ret;
    QFile file(name);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return ret;
    bool first=true;
    int num_samples=0;
    QList<Analysis*> anas;
    int n=0;
    while(!file.atEnd()) {

        QByteArray line = file.readLine();
        QStringList splitted = QString(line).split(";");
        if(first) {
            num_samples = splitted.size()-9;
            qDebug() << "num_samples: " << num_samples;
            for(int i=9; i<splitted.size(); i++) {
                anas.append(new Analysis(splitted.at(i).simplified()));
            }
            first=false;
        } else {
            for(int i=9; i<splitted.size(); i++) {
                anas.at(i-9)->buckets[n] = sqrt(splitted.at(i).toDouble());
                //qDebug() << splitted.at(i).toDouble();

            }
            n++;
        }

        //qDebug() << "Num: " << splitted.size();
        /*QString name = splitted.at(0);
        name = name.replace("\"","");
        Analysis ana(name);
        for(int i=1; i<splitted.size(); i++) {
            ana.buckets[i-1] = splitted.at(i).toDouble();
        }
        ret.append(ana);*/
    }
    for(int i=0; i<anas.size(); i++)
        ret.append(*anas.at(i));
    qDeleteAll(anas);
    qDebug() << "number of buckets/peaks: " << n;
    return ret;
}

double calc_cosine(Analysis a, Analysis b, float thresh=0) {
    if(a.buckets.size()!=b.buckets.size()) {
        qDebug() << "Bucket size does not match!";
        return -1;
    }
    double s1=0;
    double s2=0;
    double s3=0;

    int bs=a.buckets.size();

    for(int i=0; i<bs; i++) {
        //qDebug() << "a: " << a.buckets[i] << ", b: " << b.buckets[i];
        if((a.buckets[i]>thresh)||(b.buckets[i]>thresh)) {

            s1+=a.buckets[i]*b.buckets[i];
            s2+=a.buckets[i]*a.buckets[i];
            s3+=b.buckets[i]*b.buckets[i];
        }
    }

    double c = s1/sqrt(s2*s3);
    //qDebug() << "s1: " << s2 << ", s2: " << s2 << ", s3: " << s3 << ", cosine: " << c << ", bs: " << bs;
    return c;
}

bool write_cosine_file(QList<Analysis> analysis, QString name, float thresh) {
    QFile outfile(name);
    if(outfile.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream ts(&outfile);
        //ts << "source\ttarget\tcosine" << endl;
        bool first=true;
        int num_ana = analysis.size();
        for(int i=0; i<(num_ana); i++) {
            if(first)
                first = false;
            else
                ts << ";";
            ts << analysis.at(i).name;
        }
        ts << endl;

        for(int i=0; i<(num_ana); i++) {
            bool first=true;
            ts << analysis.at(i).name << ";";
            for(int n=0; n<num_ana; n++) {
                double cos = calc_cosine(analysis.at(i),analysis.at(n),thresh);
                //if(cos>thresh) {

                    if(first)
                        first = false;
                    else
                        ts << ";";
                    ts << cos;
                //}
            }
            ts << endl;
            //qDebug() << "Finished node " << i;
        }
    } else
        return false;
    return true;
}
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    if(argc<2) {
        qInfo() << "Usage: " << argv[0] << " <inputTable> <absolute threshold (optional)>";
        return -1;
    }
    double thresh=0;
    if(argc>2)
        thresh = QString(argv[2]).toDouble();


    //QString infile = argv[1];
    //double thresh=0;
    //QString infile = "F://A-Crude-Extracts-pos//Dow_01//Heatmap_BDAL//2017-02-01_Dow_01_50000_MaxIso-b.csv";
    QFileInfo fi(infile);
    QString outfile = fi.absoluteDir().absoluteFilePath("peakListCosine.csv");
    //QDir dir(infile);
    qInfo() << "Reading file " << infile << "...";
    QList<Analysis> analysis = read_file_xcms(infile);
	//QList<Analysis> analysis = read_file_profile_analysis(infile);
    //QList<Analysis> analysis = read_file_metaboscape(infile);
    qInfo() << "Writing file " << outfile << "...";
    write_cosine_file(analysis,outfile,thresh);
    qInfo() << "Finished";

    return 0;
}