#!/usr/bin/env Rscript

# Send an email notification when jobs are done.
# Uses mailR to send, and keyringr to decrypt password. Instructions on how to
# use keyringr are found at <https://bit.ly/2AHUj0A>.
# From Moon & Stubbs <Title>, <doi>

library(keyringr)
library(mailR)

# Contents:
#   - get_os
#   - sendEmailNotification
#   - sendEmailNotificationWin
#   - sendEmailNotificationMac
#   - sendEmailNotificationLinux

# get credential
credential_label <- "UoB_outlook"

# email to use
email_add <- "glbcm@bristol.ac.uk"
host_name <- "smtp.office365.com"

get_os <- function() {
  # Gets the host OS. Script taken from <https://bit.ly/2AHUj0A>.
  #
  # Returns:
  #   The host OS type: "windows", "osx", or "linux".
  # get from Sys.info()
  sysinf <- Sys.info()

  # return value from Sys.info
  if (!is.null(sysinf)){
    os <- sysinf['sysname']

    # replace "Darwin" with "osx"
    if (os == 'Darwin') os <- "macos"
  } else {
    # find from .Platform if Sys.info returns NULL
    os <- .Platform$OS.type

    # replace "Darwin" with "osx"
    if (grepl("^darwin", R.version$os)) os <- "macos"

    # replace "linux-gpu" with "linux"
    if (grepl("linux-gnu", R.version$os)) os <- "linux"
  }

  # make result lowercase
  tolower(os)
}

sendEmailNotification <- function (subject = "Notification",
                                   body = "Job complete") {
  # Gets the host OS then chooses the correct function to send an email.
  #
  # Args:
  #   subject: text for the email subject field, in UTF-8.
  #   body:    content for the email body, in UTF-8.
  #
  # Returns:
  #   An email is sent to glbcm@bristol.ac.uk. Also prints the Java-Object key
  #   to the terminal on success.
  # choose function based on OS
  if (get_os() == "windows") {
    # send from Windows
    sendEmailNotificationWin(subject = subject,
                             body = body)
  } else if (get_os() == "macos") {
    # send from macOS
    sendEmailNotificationMac(subject = subject,
                             body = body)
  } else if (get_os() == "linux") {
    # send from linux
    sendEmailNotificationLinux(subject = subject,
                               body = body)
  }
}

sendEmailNotificationWin <- function (subject, body) {
  # Send a completed-job email from and to email_add address on Windows. The
  # password should be encrypted using Microsoft Data Protection API, and will
  # be decrypted here on-the-fly. See instructions at <https://bit.ly/2AHUj0A>.
  #
  # Args:
  #   subject: text for the email subject field, in UTF-8.
  #   body:    content for the email body, in UTF-8.
  #
  # Returns:
  #   An email is sent to glbcm@bristol.ac.uk. Also prints the Java-Object key
  #   to the terminal on success.
  # get credential path
  # path for Windows password location
  credential_path <- paste0(Sys.getenv("USERPROFILE"),
                            "\\DPAPI\\passwords\\",
                            Sys.info()["nodename"],
                            "\\",
                            credential_label,
                            ".txt")

  # send email via Outlook
  send.mail(from = email_add,
            to = email_add,
            subject = subject,
            body = body,
            encoding = "utf-8",
            authenticate = TRUE,
            smtp = list(host.name = host_name,
                        port = 587,
                        user.name = email_add,
                        passwd = decrypt_dpapi_pw(credential_path),
                        tls = TRUE))
}

sendEmailNotificationMac <- function (subject, body) {
  # Send a completed-job email from and to email_add address on macOS. The
  # password should be held in the macOS Keychain. See instructions at
  # <https://bit.ly/2AHUj0A>.
  #
  # Args:
  #   subject: text for the email subject field, in UTF-8.
  #   body:    content for the email body, in UTF-8.
  #
  # Returns:
  #   An email is sent to glbcm@bristol.ac.uk. Also prints the Java-Object key
  #   to the terminal on success.
  # send email via Outlook
  send.mail(from         = email_add,
            to           = email_add,
            subject      = subject,
            body         = body,
            encoding     = "utf-8",
            authenticate = TRUE,
            smtp         = list(host.name = host_name,
                                port      = 587,
                                user.name = email_add,
                                passwd    = decrypt_kc_pw(credential_label),
                                tls       = TRUE))
}

sendEmailNotificationLinux <- function (subject, body) {
  # Send a completed-job email from and to email_add address on Linux. The
  # password should be held in the Gnome Keyring. See instructions at
  # <https://bit.ly/2AHUj0A>.
  #
  # Args:
  #   subject: text for the email subject field, in UTF-8.
  #   body:    content for the email body, in UTF-8.
  #
  # Returns:
  #   An email is sent to glbcm@bristol.ac.uk. Also prints the Java-Object key
  #   to the terminal on success.
  # send email via Outlook
  send.mail(from = email_add,
            to = email_add,
            subject = subject,
            body = body,
            encoding = "utf-8",
            authenticate = TRUE,
            smtp = list(host.name = host_name,
                        port = 587,
                        user.name = email_add,
                        passwd = decrypt_gk_pw(credential_label),
                        tls = TRUE))
}