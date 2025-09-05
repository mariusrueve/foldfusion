%%
%% This is file `master.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% master.dtx  (with options: `class')
%% 
%% This is a generated file.
%% 
%% Copyright (C) 2010 by Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
%% 
%% This file may be distributed and/or modified under the conditions of
%% the LaTeX Project Public License, either version 1.2 of this license
%% or (at your option) any later version.  The latest version of this
%% license is in:
%% 
%%    http://www.latex-project.org/lppl.txt
%% 
%% and version 1.2 or later is part of all distributions of LaTeX version
%% 1999/12/01 or later.
%% 
\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{master}
    [2016/03/10 v1.2 .dtx master thesis class]

\newif\ifmaster@oldopts\master@oldoptsfalse
\DeclareOption{old}{\master@oldoptstrue}
\newif\ifmaster@ger\master@gerfalse
\newif\ifmaster@bsc\master@bscfalse
\DeclareOption{ger}{\master@gertrue}
\DeclareOption{bsc}{\master@bsctrue}
\DeclareOption*{%
  \PassOptionsToClass{\CurrentOption}{scrreprt}%
}
\ProcessOptions\relax

\ifmaster@oldopts
  \PassOptionsToClass{bibtotoc,%
                      liststotoc}{scrreprt}%
\else
  \PassOptionsToClass{bibliography=totoc,%
                      listof=totoc}{scrreprt}%
\fi
\LoadClass[12pt,%
           openright,%
           cleardoublepage=empty,%
           parskip=half-,%
           parindent=0pt,%
           BCOR=10mm,%
           DIV=calc]{scrreprt}
\RequirePackage{appendix}
\RequirePackage[utf8]{inputenc}
\ifmaster@ger
  \RequirePackage[english,ngerman]{babel}
\else
  \RequirePackage[ngerman,english]{babel}
\fi
\RequirePackage{graphicx}
\RequirePackage{xcolor}
\RequirePackage{xspace}
\RequirePackage[automark,footsepline,plainfootsepline,headsepline]{scrlayer-scrpage}
\RequirePackage{ntheorem}
\RequirePackage[small]{caption}
\RequirePackage[sc]{mathpazo}
\linespread{1.05}
\KOMAoptions{DIV=last}
\RequirePackage[pdftex,
                bookmarks,
                bookmarksopen,
                bookmarksopenlevel=1,
                pdfstartview=FitH,
                plainpages=false,
                pdfpagelabels,
                pdfborder={0 0 0.25}]
                {hyperref}
\RequirePackage{csquotes}
\RequirePackage[backend=biber,style=numeric,sorting=none,natbib=true]{biblatex}
\RequirePackage{array}
\ifmaster@ger
  \ifmaster@bsc
    \subject{Bachelorarbeit}
  \else
    \subject{Masterarbeit}
  \fi
\else
  \ifmaster@bsc
    \subject{Bachelor Thesis}
  \else
    \subject{Master Thesis}
  \fi
\fi
\ifmaster@ger
  \newcommand{\ZBH}{
    Zentrum f{\"u}r Bioinformatik (ZBH)\\
    Universit{\"a}t Hamburg\\
    Hamburg, Deutschland\xspace}
\else
  \newcommand{\ZBH}{
    Center for Bioinformatics (ZBH)\\
    University of Hamburg\\
    Hamburg, Germany\xspace}
\fi

\newcommand{\AdvisorA}[1]{\gdef\@AdvisorA{#1}}
\newcommand{\AdvisorB}[1]{\gdef\@AdvisorB{#1}}
\newcommand{\MatrikelNo}[1]{
  \gdef\@Matrikel{#1}}
\newcommand{\AuthorAdress}[2]{
  \gdef\@AuthorStreet{#1}
  \gdef\@AuthorZip{#2}}
\newcommand{\Degree}[1]{
  \gdef\@Degree{#1}}

\ifmaster@bsc
  \Degree{B.Sc.\ Computing in Science}
\else
  \ifmaster@ger
    \Degree{Master of Science Bioinformatik}
  \else
    \Degree{Master of Science Bioinformatics}
  \fi
\fi
\newcommand{\Maketitle}{%
  \begin{titlepage}
    \begin{center}
      \vspace{1.5cm}
      \includegraphics[height=2cm]{./figures/UHH.pdf}
      \hfill
      \hspace{22.98mm}
      \\[2cm]
      {\Large {\color{red}{\@subject}}}\\[0.9cm]
      \linespread{0.8}
      {\Large\vspace{0.3cm}\@title\vspace{0.3cm}\par}
      \linespread{1.05}
      \vspace{\stretch{1}}
\ifmaster@ger
      \begin{large}
      \begin{tabular}{c}
            Vorgelegt von\\[3mm]
                \@author\\
                Matrikelnummer \@Matrikel\\[8mm]
                %\@AuthorStreet\\
                %\@AuthorZip\\[8mm]
            im Studiengang\\[3mm]
            \@Degree\\
            am Zentrum f{\"u}r Bioinformatik Hamburg\\
            und Fachbereich Informatik\\
            der MIN-Fakult{\"a}t\\[8mm]
            \@date
      \end{tabular}
      \end{large}
\else
      \begin{large}
      \begin{tabular}{c}
            Submitted by\\[3mm]
              \@author\\
              Student number \@Matrikel\\[8mm]
              %\@AuthorStreet\\
              %\@AuthorZip\\[8mm]
            in the study program\\[3mm]
            \@Degree\\
            at the Center for Bioinformatics Hamburg\\
            and the Department of Informatics\\
            of the MIN Faculty\\[8mm]
            \@date
      \end{tabular}
      \end{large}
\fi
      \vspace{\stretch{1}}
      \vfill
      \begin{large}
\ifmaster@ger
      \begin{tabular}{r@{~}l}
                 1.\ Pr{\"u}fer/in:&\@AdvisorA\\
                 2.\ Pr{\"u}fer/in:&\@AdvisorB
      \end{tabular}
\else
      \begin{tabular}{r@{~}l}
                 1.\ Examiner:&\@AdvisorA\\
                 2.\ Examiner:&\@AdvisorB
      \end{tabular}
\fi
      \end{large}
    \end{center}
  \end{titlepage}}
\newcommand{\Assertion}[0]{
  \selectlanguage{ngerman}
  \chapter*{Eidesstattliche Versicherung}
  \pagestyle{empty}
  \vspace{1\baselineskip}
Hiermit versichere ich an Eides statt, dass ich vorliegende
\ifmaster@bsc
  Bachelorarbeit im Studiengang Computing in Science
\else
  Masterarbeit im Studiengang Bioinformatik
\fi
selbstst{\"a}ndig verfasst und keine anderen als die angegebenen Hilfsmittel --
insbesondere keine im Quellenverzeichnis nicht benannten Internet-Quellen --
benutzt habe. Alle Stellen, die w{\"o}rtlich oder sinngem{\"a}{\ss} aus
Ver{\"o}ffentlichungen
entnommen wurden, sind als solche kenntlich gemacht. Ich versichere weiterhin,
dass ich die Arbeit vorher nicht in einem anderen Pr{\"u}fungsverfahren
eingereicht habe und die eingereichte schriftliche Fassung der auf dem
elektronischen Speichermedium entspricht.

  \vspace{1\baselineskip}

  Hamburg, den \@date

  \vspace{2\baselineskip}
  \rule{5cm}{0.4pt}\\
  \@author

  \vspace{2.5\baselineskip}
Ich bin mit einer Einstellung der
\ifmaster@bsc
  Bachelorarbeit
\else
  Masterarbeit
\fi
in den Bestand der Bibliothek
des Departments Informatik einverstanden.

  \vspace{1\baselineskip}

  Hamburg, den \@date

  \vspace{2\baselineskip}
  \rule{5cm}{0.4pt}\\
  \@author
  }
  \colorlet{chapter}{black!50}
  \addtokomafont{chapter}{\color{chapter}}
  \renewcommand*{\chapterformat}{%
    \begingroup
      \setlength{\unitlength}{1mm}%
      \begin{picture}(20,40)(0,5)%
      \setlength{\fboxsep}{0pt}%
      \put(20,15){\line(1,0){\dimexpr
      \textwidth-20\unitlength\relax\@gobble}}%
      \put(0,0){\makebox(15,20)[r]{%
      \fontsize{29\unitlength}{29\unitlength}\selectfont\thechapter
      \kern-.04em
      }}%
      \put(20,15){\makebox(\dimexpr
      \textwidth-20\unitlength\relax\@gobble,\ht\strutbox\@gobble)[l]{%
      \ \normalsize\color{black}\chapapp~\thechapter\autodot
      }}%
      \end{picture} %
    \endgroup
  }
  \newcommand{\Chapter}[1]{
    \chapter{#1}
    \vspace{1cm}}
  \newcommand{\ChapterNoToc}[1]{
    \chapter*{#1}
    \vspace{1cm}}
\addtokomafont{chapter}{\normalfont}
\addtokomafont{section}{\normalfont}
\addtokomafont{subsection}{\normalfont}
\addtokomafont{title}{\normalfont}
\addtokomafont{sectioning}{\normalfont}
\setkomafont{sectioning}{\normalfont\normalcolor\bfseries}
\setkomafont{descriptionlabel}{\normalfont\normalcolor\bfseries}
\theorembodyfont{\rmfamily}
\theoremstyle{definition}
\newtheorem{Theorem}{Theorem}
\newtheorem{Satz}{Satz}
\newtheorem{Lemma}{Lemma}
\newtheorem{Beispiel}{Beispiel}
\newtheorem{Definition}{Definition}

\renewcommand{\topfraction}{0.9}
\renewcommand{\bottomfraction}{0.8}
\setcounter{topnumber}{2}
\setcounter{bottomnumber}{2}
\setcounter{totalnumber}{2}
\setcounter{dbltopnumber}{2}
\renewcommand{\dbltopfraction}{0.9}
\renewcommand{\textfraction}{0.07}
\renewcommand{\floatpagefraction}{0.7}
\renewcommand{\dblfloatpagefraction}{0.7}
\pagestyle{scrheadings}
\automark[chapter]{chapter}
\automark*[section]{}
\endinput
%%
%% End of file `master.cls'.
