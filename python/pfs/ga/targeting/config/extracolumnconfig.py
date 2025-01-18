from .config import Config, Lambda

class ExtraColumnConfig(Config):
    def __init__(self,
                 lambda_func: Lambda = None):
        
        # Constant value to generate the column
        self.constant = None

        # Lambda function to generate the column
        self.lambda_func = lambda_func

        # Columns to be passed to the lambda function
        self.lambda_args = None

        # String format pattern to use to generate the column
        self.pattern = None

        # Data type of the column
        self.dtype = None
