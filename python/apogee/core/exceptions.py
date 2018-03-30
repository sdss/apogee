# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import print_function, division, absolute_import


class ApogeeError(Exception):
    """A custom core Apogee exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(ApogeeError, self).__init__(message)


class ApogeeNotImplemented(ApogeeError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(ApogeeNotImplemented, self).__init__(message)


class ApogeeApiError(ApogeeError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Apogee API'
        else:
            message = 'Http response error from Apogee API. {0}'.format(message)

        super(ApogeeAPIError, self).__init__(message)


class ApogeeApiAuthError(ApogeeAPIError):
    """A custom exception for API authentication errors"""
    pass


class ApogeeMissingDependency(ApogeeError):
    """A custom exception for missing dependencies."""
    pass


class ApogeeWarning(Warning):
    """Base warning for Apogee."""
    pass


class ApogeeUserWarning(UserWarning, ApogeeWarning):
    """The primary warning class."""
    pass


class ApogeeSkippedTestWarning(ApogeeUserWarning):
    """A warning for when a test is skipped."""
    pass


class ApogeeDeprecationWarning(ApogeeUserWarning):
    """A warning for deprecated features."""
    pass

