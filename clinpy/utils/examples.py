# simple filter of multiple conditions can replace and with or
filters = {"and": [
    ("tablename", 'field_name', '==', 'field_value'),
    ('tablename', 'field_2_name', '!=', 'field_2_value'),
]}

# nested filter with multiple conditions

filters2 = {"and": [
    ("tablename", 'field_name', '==', 'field_value'),
    {"or": [
        ("tablename", 'field_3_name', 'in', 'field_value'),
        ('tablename', 'field_4_name', 'ilike', 'field_2_value'),
    ]},
    ('tablename', 'field_2_name', '!=', 'field_2_value'),

]}

