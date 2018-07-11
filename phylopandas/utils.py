import random
import string

def get_random_id(length):
    """Generate a random, alpha-numerical id."""
    alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
    return ''.join(random.choice(alphabet) for _ in range(length))
