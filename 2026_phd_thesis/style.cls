\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{style}[2026/09/22 PhD thesis class]

% ---------------------------------------------------------------- Options ----
\newif\if@print \@printfalse
\DeclareOption{print}{\@printtrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions\relax
\LoadClass{book}

% --------------------------------------------------------------- Packages ----
\RequirePackage{setspace}
\RequirePackage{graphicx}
\graphicspath{{Figs/Raster/}{Figs/}}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{amsmath,amssymb}

% ------------------------------------------------------------- Page layout ---
\if@print
  \RequirePackage[a4paper,hmarginratio=1:1,vmarginratio=1:1,
                  scale=0.75,bindingoffset=5mm]{geometry}
\else
  \RequirePackage[a4paper,hmarginratio=1:1,vmarginratio=1:1,
                  scale=0.75]{geometry}
\fi

\renewcommand\baselinestretch{1.2}
\clubpenalty=10000
\widowpenalty=10000

% Blank pages inserted by \cleardoublepage get no header
\renewcommand{\cleardoublepage}{\clearpage\if@twoside\ifodd\c@page\else
  \hbox{}\thispagestyle{empty}\newpage
  \if@twocolumn\hbox{}\newpage\fi\fi\fi}

% ----------------------------------------------------------- Chapter style ---
\RequirePackage[Lenny]{fncychap}
\ChNameVar{\fontsize{18}{20}\usefont{OT1}{ptm}{m}{n}\selectfont}
\ChTitleVar{\huge\bfseries\fontfamily{ptm}\selectfont}

% ------------------------------------------------------------------ Headers --
\RequirePackage{fancyhdr}
\setlength{\headheight}{14.5pt}
\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{\markboth{#1}{}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}
\fancyhf{}
\fancyhead[LO]{\nouppercase\rightmark}
\fancyhead[LE,RO]{\bfseries\thepage}
\fancyhead[RE]{\nouppercase\leftmark}
\fancypagestyle{plain}{\fancyhead{}\renewcommand{\headrulewidth}{0pt}}

% ------------------------------------------------------------ Referencing ----
\RequirePackage[giveninits=true,maxbibnames=10,url=false,isbn=false]{biblatex}
\DeclareFieldFormat*{titlecase}{\MakeSentenceCase{#1}}
\DeclareNameAlias{author}{family-given}
\renewbibmacro{in:}{}

% --------------------------------------------------------------- Hyperref ----
\RequirePackage[unicode=true]{hyperref}
\RequirePackage{url}
\if@print
  \def\PHD@linkcolor{black}
\else
  \def\PHD@linkcolor{blue}
\fi
\hypersetup{
  colorlinks=true,
  linkcolor=\PHD@linkcolor, urlcolor=\PHD@linkcolor,
  citecolor=\PHD@linkcolor, anchorcolor=\PHD@linkcolor,
  bookmarksopen=true, bookmarksnumbered=true,
  breaklinks=true, linktocpage, plainpages=false,
}
\AtBeginDocument{\hypersetup{pdftitle={\@title}, pdfauthor={\@author}}}

% ------------------------------------------------------ Table of contents ----
\renewcommand{\contentsname}{Table of Contents}
\let\PHD@tableofcontents\tableofcontents
\renewcommand{\tableofcontents}{%
  \cleardoublepage
  \pdfbookmark[0]{\contentsname}{Contents}%
  \PHD@tableofcontents
  \addtocontents{toc}{\protect\vspace{-4em}}}

% ------------------------------------------------------ Title page fields ----
\newcommand{\@subtitle}{}
\newcommand{\subtitle}[1]{\renewcommand{\@subtitle}{#1}}
\newcommand{\@degreedate}{\the\year}
\newcommand{\degreedate}[1]{\renewcommand{\@degreedate}{#1}}
\newcommand{\@subjectarea}{}
\newcommand{\subjectarea}[1]{\renewcommand{\@subjectarea}{#1}}
\newcommand{\@degreetitle}{}
\newcommand{\degreetitle}[1]{\renewcommand{\@degreetitle}{#1}}
\newcommand{\@dept}{}
\newcommand{\dept}[1]{\renewcommand{\@dept}{#1}}
\newcommand{\@university}{}
\newcommand{\university}[1]{\renewcommand{\@university}{#1}}
\newcommand{\@crest}{}
\newcommand{\crest}[1]{\renewcommand{\@crest}{#1}}
\newcommand{\@supervisor}{}
\newcommand{\supervisor}[1]{\renewcommand{\@supervisor}{#1}}
\newcommand{\@supervisorrole}{Supervisors:}
\newcommand{\supervisorrole}[1]{\renewcommand{\@supervisorrole}{#1}}
\newcommand{\@candidaterole}{Candidate:}
\newcommand{\candidaterole}[1]{\renewcommand{\@candidaterole}{#1}}

% ------------------------------------------------------------- Title page ----
\renewcommand{\maketitle}{%
  \pagenumbering{alph}% unique page anchors before \frontmatter
  \thispagestyle{empty}
  \begin{singlespace}
  \centering
  % Title between two rules
  {\LARGE\bfseries
    \rule{0.98\textwidth}{0.4pt}\par\vspace{0.5em}
    \@title\par
    \ifx\@subtitle\@empty\else{\Large\mdseries\@subtitle\par}\fi
    \vspace{0.5em}\rule{0.98\textwidth}{0.4pt}\par}
  \vspace*{\stretch{.25}}
  {\@crest\par}
  \vspace*{\stretch{.15}}
  \parbox{0.5\textwidth}{\centering
    {\Large\@candidaterole\par}{\large\@author\par}}\par
  \vspace{1.5em}
  \ifx\@supervisor\@empty\else
    \parbox{0.5\textwidth}{\centering
      {\Large\@supervisorrole\par}{\large\@supervisor\par}}\par
    \vspace{1.5em}
  \fi
  {\large\@dept\par\@university\par}
  \vspace*{\stretch{.2}}
  {\small\scshape
    A thesis submitted in fulfilment of the requirements for the degree of\par
    \@degreetitle\ in \@subjectarea,\par
    \@university, \@degreedate.\par
    This thesis is for examination purposes only and is confidential to
    the examination process.\par}
  \vspace*{\stretch{.4}}
  \end{singlespace}
  \newpage}

% ---------------------------------------------- Abstract / Acknowledgements --
\newenvironment{abstract}{%
  \cleardoublepage
  \chapter*{\centering\Large\vspace{-5em}Abstract\vspace{-2em}}%
  \thispagestyle{empty}}{}

\newenvironment{acknowledgements}{%
  \cleardoublepage
  \chapter*{\centering\Large\vspace{-5em}Acknowledgements\vspace{-2em}}%
  \thispagestyle{empty}}{}

\endinput
