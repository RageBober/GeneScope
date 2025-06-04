from sklearn.feature_extraction.text import CountVectorizer

def encode_sequences(sequences, k=3):
    """
    Кодирование последовательностей ДНК/РНК в числовую форму с использованием k-mers.
    """
    vectorizer = CountVectorizer(analyzer='char', ngram_range=(k, k))
    encoded = vectorizer.fit_transform(sequences).toarray()
    feature_names = vectorizer.get_feature_names_out()
    return encoded, feature_names
