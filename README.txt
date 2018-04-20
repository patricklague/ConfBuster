============================================================
                   ConfBuster Web Server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Servers Installation and Deployment Instructions
============================================================


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                           About
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
File Author : Gabriel Bégin
Contact     : http://confbuster.org

File        : README.txt
Encoding    : UTF-8


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                       Please Cite !
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bégin, G., Barbeau, X., Vincent, A.T. & Lagüe, P. 
ConfBuster Web Server: a free web application for macrocycle
conformational search and analysis (in preparation).


Barbeau, X., Vincent, A.T. & Lagüe, P., (2018).
ConfBuster: Open-Source Tools for Macrocycle Conformational 
Search and Analysis. 
Journal of Open Research Software. 6(1),p.1. 
DOI: http://doi.org/10.5334/jors.189


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     Table of Contents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(01-A) Please read before installing - Architecture
(01-B) Please read before installing - System Commands

(02-A) External Dependencies - ConfBuster
(02-B) External Dependencies - Python
(02-C) External Dependencies - Apache
(02-D) External Dependencies - PHP
(02-E) External Dependencies - MySQL MariaDB

(03-A) Web Server - Installation
(03-B) Web Server - Deployment
(03-C) Web Server - Configuration

(04-A) Job Preparation Server - Installation
(04-B) Job Preparation Server - Deployment
(04-C) Job Preparation Server - Configuration

(05-A) Compute Server - Installation
(05-B) Compute Server - Deployment
(05-C) Compute Server - Configuration


------------------------------------------------------------
(01-A) Please read before installing - System Architecture
------------------------------------------------------------
The ConfBuster Web Server system was made to be 
installed on two different physical computers :

Computer 1 : * WebServer

Computer 2 : * JobPreparation Server
             * Compute Server

However, it is possible to install the three software
servers on the same computer.

Please refer to the UML Architecture Diagram located in the 
ConfBuster Web Server distribution for additional details.


------------------------------------------------------------
(01-B) Please read before installing - System Commands
------------------------------------------------------------
This installation/deployment protocol was tested on :

Linux Mint 18.3 Cinnamon 3.6.7 64-bits



------------------------------------------------------------
(02-A) External Dependencies - ConfBuster
------------------------------------------------------------
* Install ConfBuster
  - The lastest version
  
  https://github.com/patricklague/ConfBuster

  At the end of the installation, ConfBuster should be in
  the system $PATH (eg. /usr/bin) with all its dependencies 
  installed.


------------------------------------------------------------
(02-B) External Dependencies - Python
------------------------------------------------------------
* Install Python 2
  - Version >= 2.7.12
  
  sudo apt install python


* Install Python 2 MySQL Module
  - Version >= 1.3.7
  
  sudo apt install python-mysqldb


------------------------------------------------------------
(02-C) External Dependencies - Apache
------------------------------------------------------------
* Install the Apache HTTP Server 
  - Version >= 2.4.18

  sudo apt install apache2


* Install the Apache PHP Module
  
  sudo apt install libapache2-mod-php7.0


------------------------------------------------------------
(02-D) External Dependencies - PHP
------------------------------------------------------------
* Install PHP-5
  - Version (PHP-5) >= 5.4.45

  sudo apt install php7.0
  

* Install PHP MySQL Module
  
  sudo apt install php7.0-mysql


------------------------------------------------------------
(02-E) External Dependencies - MySQL MariaDB
------------------------------------------------------------
* Install MySQL MariaDB Database Server
  - MySQL Version >= 15.1
  - MariaDB Version >= 10.0.33

  sudo apt install mariadb-server-10.0


------------------------------------------------------------
(03-A) Web Server - Installation
------------------------------------------------------------
..i) Create a folder named "ConfBuster" into the "/var/www" 
     directory.

     sudo mkdir /var/www/ConfBuster


.ii) Move the "WebServer" folder of the GitHub Distribution 
     into the directory "/var/www/ConfBuster".

     sudo mv WebServer /var/www/ConfBuster/


------------------------------------------------------------
(03-B) Web Server - Deployment
------------------------------------------------------------
...i) Open the Apache configuration file.

      sudo nano /etc/apache2/apache2.conf


..ii) Configure the home page of the web site by adding the
      following lines at the end of the file:

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #                Directory Index
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # List of all possible file name that can be 
      # used as a home page.
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      <IfModule dir_module>
          DirectoryIndex index.html index.php \
                         Home.html home.html home.php
      </IfModule>


.iii) Open the Apache port configuration file.

      sudo nano /etc/apache2/ports.conf


..iv) Configure the port if needed (default 80).

      Listen 80


...v) Open the Apache Virtual Server configuration file.

      sudo nano /etc/apache2/sites-available/000-default.conf


..vi) Configure the port used by the Virtual server match the
      one configured in step iv) by making sure the 
      "VirtualHost" start tag as the following syntax:

      <VirtualHost *:80>


.vii) Configure the folder path of the web site by changing 
      the "DocumentRoot" statement by this one :

      DocumentRoot /var/www/ConfBuster/WebServer


viii) Reboot the Apache server

      sudo /etc/init.d/apache2 restart


------------------------------------------------------------
(03-C) Web Server - Configuration
------------------------------------------------------------
.i) Open the script used for redirecting the user's data to
    the Job Preparation Server :

    sudo nano /var/www/ConfBuster/WebServer/DataReception.php


ii) Edit the URL corresponding to the location of the 
    "ComputeJobPreparation.php" file in the JobPreparationServer.


    If the JobPreparationServer is on another computer :
 
    http://***.***.***.***:8080/BusinessLayer/ComputeJobPreparation.php


    If the JobPreparationServer is on the same computer :

    http://127.0.0.1:8080/BusinessLayer/ComputeJobPreparation.php


------------------------------------------------------------
(04-A) Job Preparation Server - Installation
------------------------------------------------------------
.i) Create a folder named "ConfBuster" into the "/var/www" 
    directory.

    sudo mkdir /var/www/ConfBuster


ii) Move the "JobPreparationServer" folder of the 
    GitHub Distribution into the directory "/var/www/ConfBuster"

    sudo mv JobPreparationServer /var/www/ConfBuster/


------------------------------------------------------------
(04-B) Job Preparation Server - Deployment
------------------------------------------------------------

..i) Open the Apache port configuration file

     sudo nano /etc/apache2/ports.conf


.ii) Add the following line at the end of the file

     Listen 8080


iii) Open the Apache Virtual Server configuration file

     sudo nano /etc/apache2/sites-available/000-default.conf


.iv) Configure the port used by the Virtual server by making
     sure the "VirtualHost" start tag as the following syntax :

     <VirtualHost *:8080>


..v) Configure the folder path of the web site by changing 
     the "DocumentRoot" statement by this one :

     DocumentRoot /var/www/ConfBuster/JobPreparationServer


.vi) Reboot the Apache server.

     sudo /etc/init.d/apache2 restart


------------------------------------------------------------
(04-C) Job Preparation Server - Configuration
------------------------------------------------------------
..i) Open the file script used for sending email :

     sudo nano /var/www/ConfBuster/JobPreparationServer/BusinessLayer/Mail/MailManager.py


.ii) Edit the Mailling information located in the 
     private attributes of the MailManager class.


iii) Open the database manager :

     nano /var/www/ConfBuster/JobPreparationServer/DatabaseLayer/DatabaseManager.php


.iv) Edit the database connection credentials located in the 
     private attributes of the DatabaseManager class.


------------------------------------------------------------
(05-A) Compute Server - Installation
------------------------------------------------------------

Since the server has to be started manually, it does not
matter where the folder is located. For convenience purpose,
it is recommended to follow this optional step :


i) Move the "ComputeServer" folder of the GitHub Distribution 
   into the directory "/var/www/ConfBuster"

   sudo mv ComputeServer /var/www/ConfBuster/


------------------------------------------------------------
(05-B) Compute Server - Deployment
------------------------------------------------------------
..i) Open the database installation script :
    
     nano /var/www/ConfBuster/ComputeServer/DatabaseLayer/Database_DDL_DML_DCL_Script.sql


.ii) Configure the credential of the user at the end of the file.


iii) Install the database with the following command :
  
     sudo mysql -u root -p < /var/www/ConfBuster/ComputeServer/DatabaseLayer/Database_DDL_DML_DCL_Script.sql


     Input upon password prompt :
  
     root


.iv) Reboot the Database server.

     sudo /etc/init.d/mysql restart


..v) Configure the Queue and the Compute Server by 
     editing the corresponding parameters located in the
     following file :

     sudo nano /var/www/ConfBuster/ComputeServer/BusinessLayer/QueueManager.py
   

.vi) Start the Compute Server.

     sudo python /var/www/ConfBuster/ComputeServer/BusinessLayer/QueueManager.py


------------------------------------------------------------
(05-C) Compute Server - Configuration
------------------------------------------------------------
..i) Open the file used for sending email :

     sudo nano /var/www/ConfBuster/ComputeServer/BusinessLayer/ManagerClasses/MailManager.py


.ii) Edit the Address of the Mail SMTP Server located in the 
     private attributes of the MailManager class.


iii) Open the database manager :

     nano /var/www/ConfBuster/ComputeServer/BusinessLayer/QueueManager.py


.iv) Edit the database connection credentials located in the 
     private attributes of the DatabaseManager class.


------------------------------------------------------------

Enjoy !

------------------------------------------------------------
End README.txt