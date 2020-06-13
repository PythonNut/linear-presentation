% ================================================================== %
%                                                                    %
%                         Requirement Things                         %
%                                                                    %
% ================================================================== %

\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{fkpaper}[2018/12/14 v0.01 paper class]

% This class based on the standard article class
\LoadClass{article}



% ================================================================== %
%                                                                    %
%                              Packages                              %
%                                                                    %
% ================================================================== %

% Sensible margins
\RequirePackage[margin=1in, tmargin=1.25in]{geometry}

% Better fonts
\RequirePackage{lmodern}

% Make sure we're inputting stuff correctly
\RequirePackage[utf8]{inputenc}

% And in english!
\RequirePackage[english]{babel}

% I've offloaded my math macros into a separate package (in this repo)
\RequirePackage{fkmath}

% For prettier code listings
\RequirePackage{fkmisc}

% For \includegraphics
\RequirePackage{graphicx}

% Fancy page header
\RequirePackage{fancyhdr}

% For multi-column environment
\RequirePackage{multicol}

% In case you want any diagrams, etc.
\RequirePackage{tikz}
\usetikzlibrary{arrows} % Better arrows
\usetikzlibrary{arrows.meta}
\usetikzlibrary{calc} % Add coordinates and stuff
\usetikzlibrary{patterns} % distinguishing shapes
\usetikzlibrary{decorations.pathreplacing} % For braces and such

% ...Which you would probably place in figure environments
\RequirePackage{float}

% For referencing page count
\RequirePackage{lastpage}

% Because I've learned my lesson about manual line breaks
\RequirePackage{parskip}

% For customizing the enumerate environment
\RequirePackage{enumitem}

% So that the whitespace around section headers can be made less
% obnoxiously large
\RequirePackage{titlesec}

% For better commands
\RequirePackage{xifthen}

% For quasi-legitimate looking tables
\RequirePackage{booktabs}

% For margin modification in subsections
\RequirePackage{changepage}

% Problem and reading set lists, for top of document
\RequirePackage{etoolbox}

% Make links clickable and colored
\RequirePackage[colorlinks=true]{hyperref}

% More intelligent cross-referencing
\RequirePackage[capitalize]{cleveref}

\NeedsTeXFormat{LaTeX2e}

% ================================================================== %
%                                                                    %
%                  Formatting for Scientific Papers!                 %
%                                                                    %
% ================================================================== %

% Use tighter margins
\newgeometry{margin=.75in, tmargin=1in}

% them appropriately sized, and what have you.
\titleformat{\section}
  {\Large\scshape\filcenter}{\scshape\thesection.}{.5em}{}
\titleformat{\subsection}
  {\large\scshape\filcenter}{\scshape\thesubsection.}{.5em}{}
\titleformat{\subsubsection}
  {\scshape\filcenter}{\bfseries\thesubsubsection.}{.5em}{}

% % Ha-ha-ha-ha, harrrr-a-ha-ha-ha-hA, ha-HA, ha-HA...
% \titlespacing{\section}{0pt}{10pt plus 4pt minus 2pt}{-5pt plus 2pt
%   minus 2pt}
% \titlespacing{\subsection}{0pt}{10pt plus 4pt minus 2pt}{-3pt plus 0pt
%   minus 0pt}
% \titlespacing{\subsubsection}{0pt}{10pt plus 4pt minus 2pt}{0pt plus
%   2pt minus 2pt}

% Make the double column format more appealing by spacing things out a
% little more
\setlength\columnsep{25pt}

% Put text lines a little closer together
\renewcommand{\baselinestretch}{.9}

% And a mild indentation
\setlength\parindent{1em}


% Change enumerate item spacing, and declare default enumerate label.
\setlist[enumerate]{itemsep=4pt, parsep=0pt, label=\arabic*)}

% block quote environment
\AtBeginEnvironment{quote}{\small}

\def\paper@title{}
\def\paper@author{}
\def\paper@author@affiliation{}
\def\paper@author@altaffiliation{}

\renewcommand{\title}[1]{\def\paper@title{#1}}
\renewcommand{\author}[1]{\def\paper@author{#1}}
\newcommand{\affiliation}[1]{\def\paper@author@affiliation{#1}}


\newcommand{\makepapertitle}{
  \thispagestyle{empty}
  \begin{center}
    \uppercase{\paper@title}

    {\scshape \small \paper@author}

    \vspace{-.3cm}
    {\footnotesize \paper@author@affiliation}
  \end{center}
}

\newcommand{\makecover}{
  \AtBeginDocument{\makepapertitle}
}

% Make headers
\newcommand{\makeheaders}{
  \pagestyle{fancy}
  \renewcommand{\headrulewidth}{0pt}
  \fancyhf{}
  \lhead{\thepage}
  \chead{\scshape \paper@author}
  \rhead{}
  \setlength{\headheight}{15.2pt}
  \setlength{\headsep}{10pt}
}


% ================================================================== %
%                                                                    %
%                          Control Sequences                         %
%                                                                    %
% ================================================================== %


% ------------------------ Creature Comforts ----------------------- %

% For referencing labeled equations w/ number in parens
\let\oldref\ref
\renewcommand{\ref}[1]{\textup{(\oldref{#1})}}

% Check if command exists
\newcommand{\checkfor}[1]{%
  \ifcsname#1\endcsname%
    Command \texttt{\textbackslash#1} is defined already. Choose
    another macro!
  \else%
    Command \texttt{\textbackslash#1} is not defined yet. Go ahead!
  \fi%
}


% ================================================================== %
%                                                                    %
%                          Misc. Formattings                         %
%                                                                    %
% ================================================================== %

% Let align environments overflow page breaks
\allowdisplaybreaks

% Define theorems and the like with a less eye-popping layout
\newtheoremstyle{snazzytheorem}  % name
    {1.5em}                      % Space above
    {0em}                       % Space below
    {\itshape}                   % Body font
    {}                           % Indent amount
    {\bfseries}                   % Theorem head font
    {.}                          % Punctuation after theorem head
    {.5em}                       % Space after theorem head
    {}  % Theorem head spec (can be left empty, meaning ‘normal’)

% Use new snazzy theorem
\theoremstyle{snazzytheorem}

% Define _____x so that we can modify qed behavior below
\newtheorem{theoremx}{Theorem}[section]
\newtheorem{corollaryx}{Corollary}[theoremx] % Reset counter w new thm
\newtheorem{lemmax}[theoremx]{Lemma} % Same counter as thm
\newtheorem{propositionx}[theoremx]{Proposition} % Same counter as thm



% End with a hollow lozenge, to visualy separate from text. Add
% \qedhere, which can be used to control where the qed gets placed,
% avoiding unnecessary line breaks and hfills.
\newenvironment{theorem}
  {\pushQED{\oldqed}\renewcommand{\qedsymbol}{$\lozenge$}\theoremx}
  {\popQED\endtheoremx}

\newenvironment{corollary}
  {\pushQED{\oldqed}\renewcommand{\qedsymbol}{$\lozenge$}\corollaryx}
  {\popQED\endcorollaryx}

\newenvironment{lemma}
  {\pushQED{\oldqed}\renewcommand{\qedsymbol}{$\lozenge$}\lemmax}
  {\popQED\endlemmax}

\newenvironment{proposition}
  {\pushQED{\oldqed}\renewcommand{\qedsymbol}{$\lozenge$}\propositionx}
  {\popQED\endpropositionx}


% And so on with definitions!
\newtheoremstyle{snazzydefinition}  % name
    {1.5em}                         % Space above
    {0em}                           % Space below
    {\rmfamily}                     % Body font
    {}                              % Indent amount
    {\bfseries}                      % Theorem head font
    {.}                             % Punctuation after theorem head
    {.5em}                          % Space after theorem head
    {}  % Theorem head spec (can be left empty, meaning ‘normal’)

\theoremstyle{snazzydefinition}

\newtheorem{definitionx}{Definition}[section]
\newtheorem{example}{Example}[section]

% Don't number remarks
\newtheorem*{remark}{Remark}

\newenvironment{definition}
  {\pushQED{\oldqed}\renewcommand{\qedsymbol}{$\lozenge$}\definitionx}
{\popQED\enddefinitionx}


% `Small' proof environment (lemmas and such)
\newenvironment{sproof}[1][Proof]
{\begin{proof}[#1]\renewcommand{\qedsymbol}{$\square$}}
{\popQED\end{proof}}



\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax



\makecover{}
\makeheaders{}

\endinput
