** Stichprobenziehung von Gemeinden mit Einwohnerzahl als Ziehungsgewicht
 * Systematische Stichprobe
 * Verschiedene VarianzschÑtzungen :
   1. g * ( y(i) - y_tot )**2   Summe g normiert auf Stichprobengrî·e
   2. Normales Ziehen mit ungleichen Wahrscheinlichkeiten mit ZurÅcklegen
** Es werden die Ergebnisse von MaxIter Stichproben in eine Datei
   geschrieben und dann ausgewertet. Insbesondere wird geprÅft, ob die
   berechneten Konfidenzintervalle korrekt Åberdecken.
*******************************************************************;
%let stich   = l.stichpr ;     * EndgÅltige Ausgabedatei             *;
%let ein     = l.gg;           * Eingabedatei                        *;
%let mos     = g_ew;           * Ziehungsgewicht, Variablenname      *;
%let y       = xxxx ;          * Variable, deren Total zu schÑtzen   *;
%let n_stp   = 250 ;           * Stichpobengrî·e                     *;
%let sum_mos = xxxxxx;         * Summe der Ziehungsgewichte          *;
%let MaxIter = 3     ;         * Anzahl Stichproben (=Iterartionen)  *;
*******************************************************************;
*------------------------------------------------------------*;
data zufall;
keep zuf;
do i=1 to &MaxIter;
  zuf = ranuni();
  output;
end;
*------------------------------------------------------------*;
%macro iterat;
%do iter=1 %to &MaxIter;
*------------------------------------------------------------*;
data Stich; * TemporÑre Ausgabdatei von einer Iteration *;
set  &ein end=Ende;
keep y_tot y_true std1 std2 kz95_1 kz95_2 kz70_1 kz70_2 Diff var1gg
  ;
label y_tot    = "Total y hochgerechnet"
      y_true   = "Total y tatsÑchlich in GG"
      dif      = "y_true - y_tot"
      var1gg   = "Varianz der GG aus Stichprobe geschÑtzt"
      std1     = "Std1: g*(y - y_tot)**2, g normiert"
      std2     = "Std2: (y*g - y_tot)**2"
      kz95_1   = "Std1*1.96 Åberdeckt y_true"
      kz95_2   = "Std2*1.96 Åberdeckt y_true"
      kz70_1   = "Std1*1.04 Åberdeckt y_true"
      kz70_2   = "Std2*1.04 Åberdeckt y_true"
   ;

array a_mos  (1,&n_stp) _temporary_; retain a_mos;
array a_var  (1,&n_stp) _temporary_; retain a_var;
array a_gew  (1,&n_stp) _temporary_; retain a_gew;

retain kum interv gezogen y_tot y_true su_gew 0 ;
*------------------------------------------------------------*;
* Ziehungsintervall und Startpunkt kum *;
if ( _n_ = 1 ) then do;
  interv = &sum_mos / &n_stp ;
  if ( interv <= 0.0001 ) then abort abend 33 ;

  set zufall point=&iter;
  kum = zuf * interv;
end; * if n=1 *;
*------------------------------------------------------------*;
* Eigentliche Ziehung *;
kum = kum + &mos;
y_true = y_true + &y ;
do while ( kum >= interv );
  kum = kum - interv ;
  gezogen = gezogen + 1;
  a_mos(gezogen) = &mos;
  a_var(gezogen) = &y;
  gewicht =  &sum_mos / &mos ; * Inverse d. Auswahlwarscheinlichkeit *;
  a_gew(gezogen) = gewicht;
  su_gew = su_gew + gewicht;

  * Hochrechnung von &y auf die Summe-&y in der Grundgesamtheit *;
  y_tot = y_tot + &y * gewicht ; ** y(i) * 1/p(i) **;
end;

if ( Ende ) then do;
  if ( gezogen ne &n_stp ) then do;
    put gezogen= ;
    abort abend 34 ;
  end;

  y_tot = y_tot / &n_stp; * weil p(i) und nicht PI(i) benutzt wurde *;
  dif  = y_true - y_tot;  * Abweichung vom wahren Wert *;

  do i=1 to &n_stp;  ** Varianz etc berechnen **;
    ** Gewichte (nur bei var1) auf Stichprobengrî·e normieren **;
    var1 = var1 +  (a_gew(i)*&n_stp/su_gew) * ( a_var(i) - y_tot )**2 ;

    * normale VarianzschÑtzung bei Ziehung mit ZurÅcklegen mit
      unterschiedlichen Auswahlwahrscheinlickkeiten *;
    var2 = var2 + ( a_gew(i)*a_var(i) - y_tot )**2 ;
  end;

  var1gg = var1 / (&n_stp - 1); *angeblich Varianz der Grundgesamtheit*;
  var1st = var1gg / &n_stp;     *angeblich Varianz des Stichprobenwerts*;
  std1   = sqrt(var1st) ;
  var2   = var2 / ( &n_stp * (&n_stp - 1) ); *Varianz des Stichp.Werts*;
  std2   = sqrt(var2) ;
end; ** if Ende **;

if ( y_tot - 1.96*std1 <= y_true <= y_tot + 1.96*std1 ) then kz95_1 = 1 ;
if ( y_tot - 1.96*std2 <= y_true <= y_tot + 1.96*std2 ) then kz95_2 = 1 ;
if ( y_tot - 1.04*std1 <= y_true <= y_tot + 1.04*std1 ) then kz70_1 = 1 ;
if ( y_tot - 1.04*std2 <= y_true <= y_tot + 1.04*std2 ) then kz70_2 = 1 ;
output stich;
run;
*------------------------------------------------------------*;
proc append;
base=&stich data=stich;
*------------------------------------------------------------*;
data _null_; run; * Macroprozessor anschubsen *;
%end; %mend; %iterat;
*******************************************************************;
*------------------------------------------------------------*;
proc print data=&stich (obs=3);
*------------------------------------------------------------*;
proc freq data=&stich;
tables kz95_1 kz95_2 kz70_1 kz70_2
  / missing;
*------------------------------------------------------------*;
proc univariate data=&stich;
var dif var1gg
  ;
*------------------------------------------------------------*;
*******************************************************************;
