from .config import Config

class ExtraColumnConfig(Config):
    def __init__(self):
        # Constant value to generate the column
        self.constant = None

        # Lambda function to generate the column
        self.lambda_func = None

        # Columns to be passed to the lambda function
        self.lambda_args = None

        # String format pattern to use to generate the column
        self.pattern = None

        # Data type of the column
        self.dtype = None
